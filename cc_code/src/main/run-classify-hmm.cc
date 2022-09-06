/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "run-classify-hmm.h"

// Standard C++ library headers:
#include <limits>
#include <string>
#include <vector>

// Local project headers:
#include "classifiers/hmm-classifier.h"
#include "common/dye-seq.h"
#include "common/radiometry.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "io/dye-seqs-io.h"
#include "io/radiometries-io.h"
#include "io/scored-classifications-io.h"
#include "main/cmd-line-out.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"
#include "util/delete.h"
#include "util/time.h"

namespace whatprot {

namespace {
using std::string;
using std::vector;
}  // namespace

void run_classify_hmm(double hmm_pruning_cutoff,
                      string dye_seqs_filename,
                      string radiometries_filename,
                      string predictions_filename) {
    double total_start_time = wall_time();

    double start_time;
    double end_time;

    start_time = wall_time();
    unsigned int num_channels;
    unsigned int total_num_dye_seqs;  // redundant, not needed.
    vector<SourcedData<DyeSeq, SourceCount<int>>> dye_seqs;
    read_dye_seqs(
            dye_seqs_filename, &num_channels, &total_num_dye_seqs, &dye_seqs);
    end_time = wall_time();
    print_read_dye_seqs(dye_seqs.size(), end_time - start_time);

    start_time = wall_time();
    unsigned int num_timesteps;
    unsigned int duplicate_num_channels;  // also get this from dye seq file.
    unsigned int total_num_radiometries;  // num radiometries across all procs.
    vector<Radiometry> radiometries;
    read_radiometries(radiometries_filename,
                      &num_timesteps,
                      &duplicate_num_channels,
                      &total_num_radiometries,
                      &radiometries);
    end_time = wall_time();
    print_read_radiometries(total_num_radiometries, end_time - start_time);

    start_time = wall_time();
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.20;
    seq_model.p_detach = 0.0055;
    seq_model.channel_models.push_back(new ChannelModel());
    seq_model.channel_models[0]->p_bleach = 0.03;
    seq_model.channel_models[0]->p_dud = 0.05;
    seq_model.channel_models[0]->bg_sig = 0.0333;
    seq_model.channel_models[0]->mu = .507;
    seq_model.channel_models[0]->sig = 0.0507;
    seq_model.channel_models[0]->stuck_dye_ratio = 0.5;  // classifiers ignore
    seq_model.channel_models[0]->p_stuck_dye_loss = 0.08;  // classifiers ignore
    seq_model.channel_models[0]->fret_eff = 0.0;
    seq_model.channel_models.push_back(new ChannelModel());
    seq_model.channel_models[1]->p_bleach = 0.06;
    seq_model.channel_models[1]->p_dud = 0.18;
    seq_model.channel_models[1]->bg_sig = 0.0333;
    seq_model.channel_models[1]->mu = .367;
    seq_model.channel_models[1]->sig = .0733;
    seq_model.channel_models[1]->stuck_dye_ratio = 0.5;  // classifiers ignore
    seq_model.channel_models[1]->p_stuck_dye_loss = 0.08;  // classifiers ignore
    seq_model.channel_models[1]->fret_eff = 0.1;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = hmm_pruning_cutoff;
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    HMMClassifier classifier(
            num_timesteps, num_channels, seq_model, seq_settings, dye_seqs);
    end_time = wall_time();
    print_built_classifier(end_time - start_time);

    start_time = wall_time();
    vector<ScoredClassification> results = classifier.classify(radiometries);
    end_time = wall_time();
    print_finished_classification(end_time - start_time);

    start_time = wall_time();
    write_scored_classifications(
            predictions_filename, total_num_radiometries, results);
    end_time = wall_time();
    print_finished_saving_results(end_time - start_time);

    double total_end_time = wall_time();
    print_total_time(total_end_time - total_start_time);
}

}  // namespace whatprot
