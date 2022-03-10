#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <vector>
#include <math.h>
#include <sstream>
#include <string>
#include <map>
#include <unistd.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <string.h>
#include <algorithm>

#define NRSAMPLES 250
#define PDF_BINS 50
#define MAX_QS 2048

#define percentile(p, n) (ceil(float(p)/100*float(n)))

struct Parameters {
    double rtt_d;
    double rtt_r;
    std::string folder;
    uint32_t n_ecn;
    uint32_t n_nonecn;
    std::string fairness;
    int nbrf;
    double link;

    Parameters() {
        rtt_d = 0;
        rtt_r = 0;
        folder = "";
        n_ecn = 0;
        n_nonecn = 0;
        fairness = "";
        nbrf = 0;
        link = 0;
    }
};

class Statistics {
  public:
    Statistics(std::string fn) {
        calculated_variance = false;
        calculated_coeffVar = false;
        _variance = NAN;
        _average = NAN;
        _coeffVar = NAN;
        _samples = NULL;
        filename_out = fn;
    }

    void samples(std::vector<double> *new_samples) {
        _samples = new_samples;
        std::sort(_samples->begin(), _samples->end());
    }

    std::vector<double> *samples() {
        return _samples;
    }

    double p(double p) {
        if (_samples != NULL && _samples->size() > 0) {
            return _samples->at(percentile(p, _samples->size()) - 1);
        }

        return NAN;
    }

    double variance() {
        if (_samples != NULL && !calculated_variance) {
            calculate_variance();
        }

        return _variance;
    }

    double average() {
        if (_samples != NULL && !calculated_variance) {
            calculate_variance();
        }

        return _average;
    }

    double coeffVar() {
        if (_samples != NULL && !calculated_coeffVar) {
            calculate_coeffVar();
        }

        return _coeffVar;
    }

    double stddev() {
        return sqrt(variance());
    }
    std::string filename_out;


  private:
    bool calculated_variance;
    bool calculated_coeffVar;
    long double _variance;
    long double _average;
    long double _coeffVar;
    std::vector<double> *_samples;

    void calculate_coeffVar() {
        if (variance() > 0 && average() > 0) {
            _coeffVar = stddev() / average();
        } else {
            _coeffVar = 0;
        }

        calculated_coeffVar = true;
    }

    void calculate_variance() {
        long double tot = 0;
        long double tot_prev = 0;
        long double sumsq = 0;
        long double sumsq_prev = 0;
        long double n_samples = (long double) _samples->size();

        for (long double val: *_samples) {
            tot_prev = tot;
            tot += val;
            if (val > 0 && tot <= tot_prev) {
                std::cout << "calc_mix: overflow in calculated_variance tot" << std::endl;
            }
            sumsq_prev = sumsq;
            sumsq += (long double) val * val;
            if (val > 0 && sumsq <= sumsq_prev) {
                std::cout << "calc_mix: overflow in calculated_variance sumsq" << std::endl;
            }
        }

        _variance = NAN;

        long double tot_sq = tot * tot;
        if (tot_sq < tot && tot > 1) {
            std::cout << "calc_mix: overflow in calculated_variance tot_sq: " << tot_sq << " tot: " << tot << std::endl;
        }
        if (n_samples > 1) {
            _variance = ((n_samples * sumsq) - (tot_sq)) / (n_samples * (n_samples - 1));
        }

        _average = tot / n_samples;
        calculated_variance = true;
    }

   
};

struct Results {
    Statistics *rate_ecn;
    Statistics *rate_nonecn;
    Statistics *win_ecn;
    Statistics *win_nonecn;
    Statistics *qs_ecn;
    Statistics *qs_nonecn;
    Statistics *drops_qs_ecn;
    Statistics *drops_qs_nonecn;
    Statistics *marks_ecn;
    Statistics *marks_nonecn;
    Statistics *util_ecn;
    Statistics *util_nonecn;
    Statistics *util; // total utilization
    std::vector<Statistics*> allstat;

    double rr_static;
    double wr_static;

    double ecn_avg;
    double nonecn_avg;

  //  uint64_t tot_sent_dropped_ecn;
    //uint64_t tot_sent_dropped_nonecn;

    Results() {
        rate_ecn = new Statistics("stat_rpf_ecn"); allstat.push_back(rate_ecn);
        rate_nonecn = new Statistics("stat_rpf_nonecn"); allstat.push_back(rate_nonecn);
        win_ecn = new Statistics("stat_win_ecn"); allstat.push_back(win_ecn);
        win_nonecn = new Statistics("stat_win_nonecn"); allstat.push_back(win_nonecn);
        qs_ecn = new Statistics("stat_qs_ecn"); allstat.push_back(qs_ecn);
        qs_nonecn = new Statistics("stat_qs_nonecn"); allstat.push_back(qs_nonecn);
        drops_qs_ecn = new Statistics("stat_drops_ecn"); allstat.push_back(drops_qs_ecn);
        drops_qs_nonecn = new Statistics("stat_drops_nonecn"); allstat.push_back(drops_qs_nonecn);
        marks_ecn = new Statistics("stat_marks_ecn"); allstat.push_back(marks_ecn);
	marks_nonecn = new Statistics("stat_marks_nonecn"); allstat.push_back(marks_nonecn);
        util_ecn = new Statistics("stat_util_ecn"); allstat.push_back(util_ecn);
        util_nonecn = new Statistics("stat_util_nonecn"); allstat.push_back(util_nonecn);
        util = new Statistics("stat_util"); allstat.push_back(util);

        rr_static = NAN;
        wr_static = NAN;

        ecn_avg = NAN;
        nonecn_avg = NAN;

    //    tot_sent_dropped_ecn = 0;
      //  tot_sent_dropped_nonecn = 0;
    }
};

struct Parameters *params = new Parameters();
struct Results *res = new Results();

void usage(int argc, char* argv[]) {
    printf("Usage: %s <folder> <e=rate_equal|d=dc_unequal> <nbr of flows per row/col> <link b/s> <rtt_d> <rtt_r> <nr ecn flows> <nr nonecn flows>\n", argv[0]);
    exit(1);
}

std::ofstream* openFileW(std::string filename) {
    std::string filename_out = params->folder + "/" + filename;
    std::ofstream *f;
    f = new std::ofstream(filename_out.c_str());

    if (!f->is_open()) {
        std::cerr << "error opening file: " << filename_out << std::endl;
        exit(1);
    }

    return f;
}

void writeToFile(Parameters *param, Statistics *s) {
    std::stringstream data;
   // data << "# num_flows average p1 p25 p75 p99 stddev" << std::endl;
    uint64_t num_flows = 0;
    if (0 == s->filename_out.compare(s->filename_out.length() - 7, 7, "_nonecn"))
        num_flows = param->n_nonecn;
    else if (0 == s->filename_out.compare(s->filename_out.length() - 4, 4, "_ecn"))
        num_flows = param->n_ecn;
    else
        num_flows = param->n_ecn + param->n_nonecn;
    data << "s" << num_flows <<  " " << s->average() << " " << s->p(1) << " " << s->p(25) << " " << s->p(75) << " " << s->p(99) << " " << s->stddev() << std::endl;
   
    std::ofstream *fs = openFileW(s->filename_out);
    *fs << data.str();
    fs->close();
}


void dmPDF(Statistics *drops_ecn, Statistics *drops_nonecn, Statistics *marks_ecn, Statistics *marks_nonecn, int i) {
    std::ofstream *f_decn_pdf = openFileW("d_pf_ecn_pdf");
    std::ofstream *f_mecn_pdf = openFileW("m_pf_ecn_pdf");
    std::ofstream *f_dnonecn_pdf = openFileW("d_pf_nonecn_pdf");
    std::ofstream *f_mnonecn_pdf = openFileW("m_pf_nonecn_pdf");

    std::vector<double> *samples_drops_ecn = drops_ecn->samples();
    std::vector<double> *samples_drops_nonecn = drops_nonecn->samples();
    std::vector<double> *samples_marks_ecn = marks_ecn->samples();
    std::vector<double> *samples_marks_nonecn = marks_nonecn->samples();

    uint32_t decn_pdf[PDF_BINS];
    uint32_t dnonecn_pdf[PDF_BINS];
    uint32_t mecn_pdf[PDF_BINS];
    uint32_t mnonecn_pdf[PDF_BINS];

    bzero(decn_pdf, sizeof(uint32_t)*PDF_BINS);
    bzero(dnonecn_pdf, sizeof(uint32_t)*PDF_BINS);
    bzero(mecn_pdf, sizeof(uint32_t)*PDF_BINS);
    bzero(mnonecn_pdf, sizeof(uint32_t)*PDF_BINS);


    uint32_t max = 0;

    if (samples_drops_ecn->back() > max)
        max = samples_drops_ecn->back();
    if (samples_drops_nonecn->back() > max)
        max = samples_drops_nonecn->back();

    uint32_t binsize = max/PDF_BINS;
    uint32_t b;

    for (double val: *samples_drops_ecn) {
        b = val / binsize;
        if (b >= PDF_BINS)
            b = PDF_BINS - 1;
        decn_pdf[b]++;
    }

    for (double val: *samples_drops_nonecn) {
        b = val / binsize;
        if (b >= PDF_BINS)
            b = PDF_BINS - 1;
        dnonecn_pdf[b]++;
    }

    if (samples_marks_ecn->back() > max)
        max = samples_marks_ecn->back();
    if (samples_marks_nonecn->back() > max)
        max = samples_marks_nonecn->back();

    binsize = max/PDF_BINS;
    for (double val: *samples_marks_ecn) {
        b = val / binsize;
        if (b >= PDF_BINS)
            b = PDF_BINS - 1;
        mecn_pdf[b]++;
    }

    for (double val: *samples_marks_nonecn) {
        b = val / binsize;
        if (b >= PDF_BINS)
            b = PDF_BINS - 1;
        mnonecn_pdf[b]++;
    }


    for (int i = 0; i < PDF_BINS; ++i) {
        *f_decn_pdf << i << " " << decn_pdf[i] << std::endl;
        *f_mecn_pdf << i << " " << mecn_pdf[i] << std::endl;
        *f_dnonecn_pdf << i << " " << dnonecn_pdf[i] << std::endl;
        *f_mnonecn_pdf << i << " " << mnonecn_pdf[i] << std::endl;
    }

    f_decn_pdf->close();
    f_mecn_pdf->close();
    f_dnonecn_pdf->close();
    f_mnonecn_pdf->close();
}

void rPDF(Statistics *rate_ecn, Statistics *rate_nonecn, char fairness, int n_ecn, int n_nonecn, int nbrf) {
    std::ofstream *f_recn_pdf = openFileW("r_pf_ecn_pdf");
    std::ofstream *f_rnonecn_pdf = openFileW("r_pf_nonecn_pdf");
    uint64_t recn_pdf[PDF_BINS];
    uint64_t rnonecn_pdf[PDF_BINS];
    bzero(recn_pdf, sizeof(uint64_t)*PDF_BINS);
    bzero(rnonecn_pdf, sizeof(uint64_t)*PDF_BINS);

    std::vector<double> *samples_rate_ecn = rate_ecn->samples();
    std::vector<double> *samples_rate_nonecn = rate_nonecn->samples();

    uint32_t max = fairness == 'e' ? n_ecn + n_nonecn : (n_ecn ? n_ecn : n_nonecn);
    if (max == 0)
        max = 1;

    max = 10000000/max/nbrf;

    uint32_t binsize = max/PDF_BINS;
    uint32_t b;

    for (double val: *samples_rate_ecn) {
        b = val / binsize;
        if (b >= PDF_BINS)
            b = PDF_BINS - 1;
        recn_pdf[b]++;
    }

    for (double val: *samples_rate_nonecn) {
        b = val / binsize;
        if (b >= PDF_BINS)
            b = PDF_BINS - 1;
        rnonecn_pdf[b]++;
    }

    for (int i = 0; i < PDF_BINS; ++i) {
        *f_recn_pdf << i << " " << recn_pdf[i] << std::endl;
        *f_rnonecn_pdf << i << " " << rnonecn_pdf[i] << std::endl;
    }

    f_recn_pdf->close();
    f_rnonecn_pdf->close();
}

void readFileMarks(std::string filename_marks, Statistics *stats, std::string filename_tot) {
    std::ifstream infile_marks(filename_marks.c_str());
    std::ifstream infile_tot(filename_tot.c_str());

    double marks;
    double tot_packets;
    std::vector<double> *samples = new std::vector<double>();

    for (int s = 0; s < NRSAMPLES; ++s) {
        for (int colnr = 0; colnr < 3; ++colnr) {
            if (infile_marks.eof() || infile_tot.eof())
                break;

            infile_marks >> marks;

            if (colnr == 2) {
                infile_tot >> tot_packets;
                double marks_perc = 0;

                if (tot_packets > 0) {
                    marks_perc = marks * 100 / tot_packets;
                    samples->push_back(marks_perc);
                }

            }
        }
    }

    infile_marks.close();
    infile_tot.close();
    stats->samples(samples);
}

void readFileDrops(std::string filename_drops, Statistics *stats, std::string filename_tot) {
    std::ifstream infile_drops(filename_drops.c_str());
    std::ifstream infile_tot(filename_tot.c_str());

    double drops;
    double tot_packets;
    std::vector<double> *samples = new std::vector<double>();

    for (int s = 0; s < NRSAMPLES; ++s) {
        for (int colnr = 0; colnr < 3; ++colnr) {
            if (infile_drops.eof() || infile_tot.eof())
                break;

            infile_drops >> drops;

            if (colnr == 2) {
                infile_tot >> tot_packets;
                double drops_perc = 0;

                if (tot_packets+drops > 0) {
                    drops_perc = drops*100/(tot_packets+drops);
                    samples->push_back(drops_perc);
                }


                if (drops_perc > 100)
                    std::cout << "too large drops perc: " << drops_perc << std::endl;
            }
        }
    }

    infile_drops.close();
    infile_tot.close();
    stats->samples(samples);
}

void readFileRate(std::string filename, int nrflows, Statistics *stats_rate, Statistics *stats_win, double avg_qs, double rtt) {
    std::ifstream infile(filename.c_str());
    double rate;

    std::vector<double> *samples_rate = new std::vector<double>();
    std::vector<double> *samples_win = new std::vector<double>();

    for (int s = 0; s < NRSAMPLES; ++s) {
        for(std::string line; getline(infile, line);) {
            std::istringstream iss(line);
            int colnr = 0;

            while (iss >> rate) {
                rate /= 8;
                if (colnr++ >= 2) {
                    double win = 0;
                    if (avg_qs != 0) {
                        win = rate*(avg_qs+rtt)/1000;
                    }

                    samples_rate->push_back(rate);
                    samples_win->push_back(win);
                }
            }

            for (int i = colnr; i < nrflows; ++i) {
                samples_rate->push_back(0);
                samples_win->push_back(0);
            }
        }
    }

    infile.close();
    stats_rate->samples(samples_rate);
    stats_win->samples(samples_win);
}

void getSamplesUtilization() {
    std::string filename_ecn = params->folder + "/rate_ecn";
    std::string filename_nonecn = params->folder + "/rate_nonecn";

    std::ifstream infile_ecn(filename_ecn.c_str());
    std::ifstream infile_nonecn(filename_nonecn.c_str());

    std::vector<double> *samples_ecn = new std::vector<double>();
    std::vector<double> *samples_nonecn = new std::vector<double>();
    std::vector<double> *samples = new std::vector<double>();
    double rate_ecn;
    double rate_nonecn;
    double util_ecn;
    double util_nonecn;
    double util;
    double link_bytes_ps = (double)params->link*125000;

   
    // each line consists of three numbers, and we only want the last number
    for (int s = 0; s < NRSAMPLES; ++s) {

        if (infile_ecn.eof() || infile_nonecn.eof()) {
            break;
        }
        for (int colnr = 0; colnr < 3; ++colnr) {
            infile_ecn >> rate_ecn;
            infile_nonecn >> rate_nonecn;
            if (colnr == 2) {
                util_ecn = rate_ecn/8 * 100 / link_bytes_ps;
                util_nonecn = rate_nonecn/8 * 100 / link_bytes_ps;
                util = (rate_ecn+rate_nonecn)/8 * 100 / link_bytes_ps;
                samples_ecn->push_back(util_ecn);
                samples_nonecn->push_back(util_nonecn);
                samples->push_back(util);
            }
        }
    }

    infile_ecn.close();
    infile_nonecn.close();
    res->util_ecn->samples(samples_ecn);
    res->util_nonecn->samples(samples_nonecn);
    res->util->samples(samples);
}

void readFileQS(std::string filename, Statistics *stats) {
    std::ifstream infile(filename.c_str());
    if (!infile.is_open()) {
        std::cerr << "calc_mix: Error opening file for reading: " << filename << std::endl;
        exit(1);
    }

    std::vector<double> *samples = new std::vector<double>();

    // Columns in file we are reading:
    // <queuing delay in us> <number of packes not dropped> <number of packets dropped>

    // we don't skip any samples for this one, as the input data
    // is already aggregated over all samples

    while (1) {
        double us;
        double nrpackets;
        double drops;

        infile >> us; /* number of us each packet represents */
        infile >> nrpackets;
        infile >> drops;

        if (infile.eof()) {
            break;
        }

        for (int i = 0; i < nrpackets; ++i) {
            samples->push_back(us/1000); // push back in ms
        }
    }

    infile.close();

    stats->samples(samples);
}

void getSamplesRateMarksDrops() {
    readFileRate(params->folder + "/flows_rate_ecn", params->n_ecn, res->rate_ecn, res->win_ecn, res->qs_ecn->average(), params->rtt_d);
    readFileMarks(params->folder + "/marks_ecn", res->marks_ecn, params->folder + "/packets_ecn");
    readFileDrops(params->folder + "/drops_ecn", res->drops_qs_ecn, params->folder + "/packets_ecn");
    readFileRate(params->folder + "/flows_rate_nonecn", params->n_nonecn, res->rate_nonecn, res->win_nonecn, res->qs_nonecn->average(), params->rtt_r);
    readFileMarks(params->folder + "/marks_nonecn", res->marks_nonecn, params->folder + "/packets_nonecn");
    readFileDrops(params->folder + "/drops_nonecn", res->drops_qs_nonecn, params->folder + "/packets_nonecn");
}

void getSamplesQS() {
    readFileQS(params->folder + "/queue_packets_drops_ecn_pdf", res->qs_ecn);
    readFileQS(params->folder + "/queue_packets_drops_nonecn_pdf", res->qs_nonecn);
}

void loadParameters(int argc, char **argv) {
    if (argc < 9) {
        usage(argc, argv);
    }

    params->folder = argv[1];
    params->fairness = argv[2];
    params->nbrf = atoi(argv[3]);
    params->link = atoi(argv[4]);
    params->rtt_d = (double) atoi(argv[5]);
    params->rtt_r = (double) atoi(argv[6]);
    params->n_ecn = atoi(argv[7]);
    params->n_nonecn = atoi(argv[8]);

    if (params->fairness.length() != 1) {
        usage(argc, argv);
    }
}

int main(int argc, char **argv) {
    loadParameters(argc, argv);

    getSamplesQS();
    getSamplesRateMarksDrops();
    getSamplesUtilization();

    if (params->n_nonecn > 0) {
        res->rr_static = res->rate_ecn->average() /  res->rate_nonecn->average();
        res->wr_static = res->win_ecn->average() / res->win_nonecn->average();
    }

    //rPDF(res->rate_ecn, res->rate_nonecn, params->fairness[0], params->n_ecn, params->n_nonecn, params->nbrf);
    //dmPDF(res->drops_qs_ecn, res->drops_qs_nonecn, res->marks_ecn, res->marks_nonecn, i);

    if (res->drops_qs_nonecn->p(99) > 100) {
        std::cerr << "too high drops p99: " << res->drops_qs_nonecn->p(99) << std::endl;
        exit(1);
    }

    for (auto it = res->allstat.begin(); it != res->allstat.end(); ++it) {
        Statistics* stat = *it;
        writeToFile(params, stat);
    }

    std::ofstream *f_avgrate_ecn = openFileW("avgrate_ecn");
    if (f_avgrate_ecn->is_open()){
        *f_avgrate_ecn << (int)res->rate_ecn->average();
        f_avgrate_ecn->close();
    }
    std::ofstream *f_avgrate_nonecn = openFileW("avgrate_nonecn");
    if (f_avgrate_nonecn->is_open()){
        *f_avgrate_nonecn << (int)res->rate_nonecn->average();
        f_avgrate_nonecn->close();
    }

    std::ofstream *f_avgqs_ecn = openFileW("avgqs_ecn");
    if (f_avgqs_ecn->is_open()){
        *f_avgqs_ecn << (int)res->qs_ecn->average();
        f_avgqs_ecn->close();
    }
    std::ofstream *f_avgqs_nonecn = openFileW("avgqs_nonecn");
    if (f_avgqs_nonecn->is_open()){
        *f_avgqs_nonecn << (int)res->qs_nonecn->average();
        f_avgqs_nonecn->close();
    }

    std::ofstream *f_rr = openFileW("rr");
 	*f_rr << res->rr_static;
 	f_rr->close();

    std::ofstream *f_wr = openFileW("wr");
    *f_wr << res->wr_static;
    f_wr->close();

    std::ofstream *f_fairwin = openFileW("fairwin");
    double fairwin = params->link*125000/(params->n_ecn/(res->qs_ecn->average()+params->rtt_d) + params->n_nonecn/(res->qs_nonecn->average()+params->rtt_r))/1000;
    *f_fairwin << fairwin;
    f_fairwin->close();



    return 0;
}
