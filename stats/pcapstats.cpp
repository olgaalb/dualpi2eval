#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <vector>
#include <math.h>
#include <sstream>
#include <string>
#include <unistd.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <string.h>
#include <algorithm>
#include <pcap.h>
#include <netinet/ip.h>
#include <netinet/if_ether.h> /* includes net/ethernet.h */
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/time.h>


typedef u_int32_t u32; // we use "kernel-style" u32 variables in numbers.h

#include "numbers.h"

#define percentile(p, n) (round(float(p)/100*float(n) + float(1)/2))
#define _PARSE(type, name, offset) type *name = (type *)(buffer + offset)
#define PARSE(type, name) _PARSE(type, name, offset)
#define QS_LIMIT 2048


void usage(int argc, char* argv[])
{
    printf("Usage: %s <pcap file> <output_folder> \n", argv[0]);
    exit(1);
}

std::string IPtoString(in_addr_t ip) {
    struct sockaddr_in ipip;
    ipip.sin_addr.s_addr = ip;
    return std::string(inet_ntoa(ipip.sin_addr));
}

 
static inline u32 qdelay_decode(u32 value)
{
    /* Input value is originally time in ns right shifted 15 times
     * to get division by 1000 and units of 32 us. The right shifting
     * by 10 to do division by 1000 actually causes a rounding
     * we correct by doing (x * (1024/1000)) here.
     */
    return fl2int(value, QDELAY_M, QDELAY_E) * 32 * 1.024;
}

double twoDecimalPlaces(double number) {
    uint32_t scaled = (uint32_t)round(number * 100);
    return ((double)number)/100;
}

int getIPHeaderOffset(pcap_t* file)
{
    int linktype = pcap_datalink(file);
    if (linktype == DLT_EN10MB) 
        return ETHER_HDR_LEN;
    else if (linktype == DLT_RAW)
        return 0;
    else if (linktype == DLT_C_HDLC)
        return 4;
    else if (linktype == DLT_NULL)
        return 4;
    else {
        std::cout << "Linktype: " << linktype << "not implemented yet." << std::endl;
        exit(1);
    }

}

struct timeval usToTs(int qdelay_us)
{
    struct timeval ts;
    bzero(&ts, sizeof(struct timeval));
    int sec = qdelay_us/1000000;
    int us = qdelay_us % 1000000;
    ts.tv_sec = sec;
    ts.tv_usec = us;
    return ts;
}

void processFile(std::string fn, std::string folder)
{    
    char errbuf[PCAP_ERRBUF_SIZE];
    const u_char *data;
    struct pcap_pkthdr h;
    in_addr_t src_ip, dst_ip;
    uint8_t tos;  
    struct iphdr* iph;
    struct timeval ts_pcap, ts_arrival, ts_delay;
    int iphdr_offset, drops, qdelay, qdelay_encoded;
    std::vector<struct timeval> lq_arrival;
    std::vector<struct timeval> lq_departure;

    int qdelay_decode_table[QS_LIMIT];


    for (int i = 0; i < QS_LIMIT; ++i) {
        qdelay_decode_table[i] = qdelay_decode(i);
    }
    pcap_t *file = pcap_open_offline(fn.c_str(), errbuf);
    if (!file)
    {
        exit(1);
    }

    iphdr_offset = getIPHeaderOffset(file);

    while (file)
    {
        data = pcap_next(file, &h);
        if (data == NULL) { 
            // write output
            return;
        }
        ts_pcap = h.ts;
       
        iph = (struct iphdr*) (data + iphdr_offset);
        if (iph->version == 6) {
            continue;
        } else {
            src_ip = iph->saddr;
            dst_ip = iph->daddr;
            uint8_t protocol = (uint8_t) iph->protocol;
            if (protocol != IPPROTO_TCP && protocol != IPPROTO_UDP)
                continue;
        }
        
        uint16_t id = ntohs(iph->id);
      
        qdelay_encoded = id & 2047; // 2047 = 0b0000011111111111
        qdelay = qdelay_decode_table[qdelay_encoded];
        
        ts_delay = usToTs(qdelay);

        timersub(&ts_pcap, &ts_delay, &ts_arrival);

        tos = ntohl(src_ip);
        if (tos & 3) {
            lq_arrival.push_back(ts_arrival);
            lq_departure.push_back(ts_pcap);
        } 
        
        if (lq_arrival.size() >= 1000)
            break;

    }

    std::sort(lq_arrival.begin(), lq_arrival.end(), [] (struct timeval a, struct timeval b) {
        return (a.tv_sec < b.tv_sec) || (a.tv_sec == b.tv_sec && a.tv_usec < b.tv_usec);
    });

    std::sort(lq_departure.begin(), lq_departure.end(), [] (struct timeval a, struct timeval b) {
        return (a.tv_sec < b.tv_sec) || (a.tv_sec == b.tv_sec && a.tv_usec < b.tv_usec);
    });

    struct timeval init = lq_arrival.front();
    struct timeval arrival_norm, arrival_norm_prev, departure_norm, departure_norm_prev;

    std::ofstream f_qs, f_qs_raw;
    std::string fn_out = folder + "/qs_l_timeline";
    std::string fn_out_raw = folder + "/qs_l_timeline_raw";

    f_qs.open(fn_out.c_str());
    if (!f_qs.is_open()) {
        std::cerr << "pcapstats: Error opening file for writing: " << fn_out << std::endl;
        exit(1);
    }

    f_qs_raw.open(fn_out_raw.c_str());
    if (!f_qs_raw.is_open()) {
        std::cerr << "pcapstats: Error opening file for writing: " << fn_out << std::endl;
        exit(1);
    }

    int qs = 0;
    auto itd = lq_departure.begin();
    auto ita = lq_arrival.begin();
    while (ita != lq_arrival.end() && itd != lq_departure.end())
    {
        timersub(ita, &init, &arrival_norm);
        timersub(itd, &init, &departure_norm);

        while (ita != lq_arrival.end() && arrival_norm.tv_usec <= departure_norm.tv_usec)
        {
            f_qs_raw << arrival_norm.tv_usec << " " << ++qs << std::endl;
            f_qs << arrival_norm.tv_usec << " " << qs << std::endl;

            arrival_norm_prev = arrival_norm;
            ita++;
            if (ita != lq_arrival.end())
                timersub(ita, &init, &arrival_norm);

            while (arrival_norm_prev.tv_usec < arrival_norm.tv_usec && arrival_norm_prev.tv_usec < departure_norm.tv_usec) {
                arrival_norm_prev.tv_usec = arrival_norm_prev.tv_usec + 1;
                f_qs << arrival_norm_prev.tv_usec << " " << qs << std::endl;
            }
        }

        while (itd != lq_departure.end() && departure_norm.tv_usec < arrival_norm.tv_usec)
        {
            f_qs_raw << departure_norm.tv_usec << " " << --qs << std::endl;
            f_qs << departure_norm.tv_usec << " " << qs << std::endl;

            departure_norm_prev = departure_norm;
            
            itd++;
            if (itd != lq_departure.end())
                timersub(itd, &init, &departure_norm);
            while (departure_norm_prev.tv_usec < departure_norm.tv_usec && departure_norm_prev.tv_usec < arrival_norm.tv_usec) {
                departure_norm_prev.tv_usec = departure_norm_prev.tv_usec + 1;
                f_qs << departure_norm_prev.tv_usec << " " << qs << std::endl;
            }

        }

        while (ita != lq_arrival.end() && departure_norm_prev.tv_usec < arrival_norm.tv_usec) {
            departure_norm_prev.tv_usec = departure_norm_prev.tv_usec + 1;
            f_qs << departure_norm_prev.tv_usec << " " << qs << std::endl;
        }

    }
    f_qs.close();
    f_qs_raw.close();

}

// transfer size/ratein bytes per second/ + 2rtt
int main(int argc, char **argv)
{
	if (argc < 3)
		usage(argc,argv);

    std::string filename = argv[1];
    std::string output_folder = argv[2];
    processFile(filename, output_folder);
    
    return 0;
}

