#include <visualization.h>
using namespace std;
namespace plt = matplotlibcpp;

// void _set_limits(map<string, float>& ylims, map<string, float>& ylims_unispg);

// void draw(vector<vector<int> >& exonIntervals, vector<vector<int> >& exonIntervals_unispg, string filename);

void _set_limits(map<string, float>& ylims, map<string, float>& ylims_unispg) {
    ylims["intron_max"] = ylims["exon_max"]*0.9;
    ylims["intron_min"] = (ylims["exon_max"] + ylims["exon_min"])/2.0;
    ylims["bar_min"] = ylims["exon_max"]+0.2;
    ylims["bar_max"] = ylims["bar_min"]+(ylims["exon_max"]-ylims["exon_min"])/5.0;

    ylims_unispg["intron_max"] = ylims_unispg["exon_max"]*0.9;
    ylims_unispg["intron_min"] = (ylims_unispg["exon_max"] + ylims_unispg["exon_min"])/2.0;
    ylims_unispg["bar_min"] = ylims_unispg["exon_max"]+0.2;
    ylims_unispg["bar_max"] = ylims_unispg["bar_min"]+(ylims_unispg["exon_max"]-ylims_unispg["exon_min"])/5.0;
}

void _transform_spans(vector<vector<int> >& exonIntervals, vector<vector<int> >& exonIntervals_unispg, float minExonLen, float minExonLen_unispg, int totalSpan, int totalSpan_unispg, int numExons, int numExons_unispg) {
    vector<float> span_lens;
    float max_len = 0.0;
    for (int i = 0; i < exonIntervals.size(); i++) {
        fprintf(stderr, "exonIntervals[%d][%d]: %d\n", i, 1,exonIntervals[i][1]);
        fprintf(stderr, "exonIntervals[%d][%d]: %d\n", i, 0,exonIntervals[i][0]);
        float span = static_cast<float>(exonIntervals[i][1] - exonIntervals[i][0]);
        fprintf(stderr, "span: %f\n", span);
        span_lens.push_back(span);
        if (span > max_len) {
            max_len = span;
        }
        fprintf(stderr, "max_len: %f\n", max_len);
    }
    fprintf(stderr, "span_lens.size(): %d\n", span_lens.size());
    vector<float> span_ratios;
    if (max_len < minExonLen) {
        fprintf(stderr, "max_len: %f\n", max_len);
        for (int i = 0; i < span_lens.size(); i++) {
            float span_ratio = span_lens[i]/max_len;
            fprintf(stderr, "span_ratio: %f\n", span_ratio);
            span_ratios.push_back(span_ratio);
        }
        float expansion_factor = static_cast<float>(totalSpan)*1e-11;
        for (int i = 0; i < 10; i++) {
            float ef = pow(2, i)*expansion_factor;
            if (max_len+ef > minExonLen) {
                expansion_factor = ef;
                break;
            }
        }

        for (int i = 0; i < exonIntervals.size(); i++) {
            float mid = static_cast<float>((exonIntervals[i][0] + exonIntervals[i][1])/2);
            float f = (expansion_factor*span_ratios[i])/2;
            if (mid+f - mid-f > minExonLen) {
                // vector<int> tmp;
                // tmp.push_back(mid-f);
                // tmp.push_back(mid+f);
                exonIntervals[i][0] = static_cast<int>(mid-f);
                exonIntervals[i][1] = static_cast<int>(mid+f);
            }
            else {
                // vector<int> tmp;
                // tmp.push_back(mid-(minExonLen/2));
                // tmp.push_back(mid+(minExonLen/2));
                // transformed_intervals.push_back(tmp);
                exonIntervals[i][0] = static_cast<int>(mid-(minExonLen/2));
                exonIntervals[i][1] = static_cast<int>(mid+(minExonLen/2));
            }
        }
    } else {
        for (int i = 0; i < numExons; i++) {
            if (span_lens[i] < minExonLen) {
                float mid = static_cast<float>((exonIntervals[i][0] + exonIntervals[i][0])/2);
                // vector<int> tmp;
                // tmp.push_back(mid-(minExonLen/2));
                // tmp.push_back(mid+(minExonLen/2));
                // transformed_intervals.push_back(tmp);
                // transformed_intervals.push_back(tmp);
                exonIntervals[i][0] = static_cast<int>(mid-(minExonLen/2));
                exonIntervals[i][1] = static_cast<int>(mid+(minExonLen/2));
            } else {
                // transformed_intervals.push_back(exonIntervals[i]);
                // exonIntervals[i][0] = mid-(minExonLen/2);
                // exonIntervals[i][1] = mid+(minExonLen/2);
            }
        }
    }

    

    
    vector<float> span_lens_unispg;
    max_len = 0.0;
    for (int i = 0; i < exonIntervals_unispg.size(); i++) {
        fprintf(stderr, "exonIntervals_unispg[%d][%d]: %d\n", i, 1,exonIntervals_unispg[i][1]);
        fprintf(stderr, "exonIntervals_unispg[%d][%d]: %d\n", i, 0,exonIntervals_unispg[i][0]);
        float span = static_cast<float>(exonIntervals_unispg[i][1] - exonIntervals_unispg[i][0]);
        fprintf(stderr, "span: %f\n", span);
        span_lens_unispg.push_back(span);
        if (span > max_len) {
            max_len = span;
        }
        fprintf(stderr, "max_len: %f\n", max_len);
    }
    vector<float> span_ratios_unispg;
    if (max_len < minExonLen_unispg) {
        fprintf(stderr, "max_len: %f\n", max_len);
        for (int i = 0; i < span_lens_unispg.size(); i++) {
            float span_ratio = span_lens_unispg[i]/max_len;
            fprintf(stderr, "span_ratio: %f\n", span_ratio);
            span_ratios_unispg.push_back(span_ratio);
        }
        float expansion_factor = static_cast<float>(totalSpan_unispg)*1e-11;
        for (int i = 0; i < 10; i++) {
            float ef = pow(2, i)*expansion_factor;
            if (max_len+ef > minExonLen_unispg) {
                expansion_factor = ef;
                break;
            }
        }

        for (int i = 0; i < exonIntervals_unispg.size(); i++) {
            float mid = static_cast<float>((exonIntervals_unispg[i][0] + exonIntervals_unispg[i][1])/2);
            float f = (expansion_factor*span_ratios_unispg[i])/2;
            if (mid+f - mid-f > minExonLen_unispg) {
                // vector<int> tmp;
                // tmp.push_back(mid-f);
                // tmp.push_back(mid+f);
                exonIntervals_unispg[i][0] = static_cast<int>(mid-f);
                exonIntervals_unispg[i][1] = static_cast<int>(mid+f);
            }
            else {
                // vector<int> tmp;
                // tmp.push_back(mid-(minExonLen_unispg/2));
                // tmp.push_back(mid+(minExonLen_unispg/2));
                // transformed_intervals.push_back(tmp);
                exonIntervals_unispg[i][0] = static_cast<int>(mid-(minExonLen_unispg/2));
                exonIntervals_unispg[i][1] = static_cast<int>(mid+(minExonLen_unispg/2));
            }
        }
    } else {
        for (int i = 0; i < numExons_unispg; i++) {
            if (span_lens_unispg[i] < minExonLen_unispg) {
                float mid = static_cast<float>((exonIntervals_unispg[i][0] + exonIntervals_unispg[i][0])/2);
                // vector<int> tmp;
                // tmp.push_back(mid-(minExonLen_unispg/2));
                // tmp.push_back(mid+(minExonLen_unispg/2));
                // transformed_intervals.push_back(tmp);
                // transformed_intervals.push_back(tmp);
                exonIntervals_unispg[i][0] = static_cast<int>(mid-(minExonLen_unispg/2));
                exonIntervals_unispg[i][1] = static_cast<int>(mid+(minExonLen_unispg/2));
            } else {
                // transformed_intervals.push_back(exonIntervals_unispg[i]);
                // exonIntervals_unispg[i][0] = mid-(minExonLen_unispg/2);
                // exonIntervals_unispg[i][1] = mid+(minExonLen_unispg/2);
            }
        }
    }
}

bool _draw_intron(map<string, float>& ylims, int span0, int span1, map<string, string>& keywords) {
    float mid = static_cast<float>((span0+span1)/2.0);
    vector<float> line_left_wrp;
    line_left_wrp.push_back(static_cast<float>(span0));
    line_left_wrp.push_back(mid);
    vector<float> ylims_left_wrp;
    ylims_left_wrp.push_back(ylims["intron_min"]);
    ylims_left_wrp.push_back(ylims["intron_max"]);

    vector<float> line_right_wrp;
    line_right_wrp.push_back(mid);
    line_right_wrp.push_back(static_cast<float>(span1));
    vector<float> ylims_right_wrp;
    ylims_right_wrp.push_back(ylims["intron_max"]);
    ylims_right_wrp.push_back(ylims["intron_min"]);

    // std::map<string, string> keywords;
    // if (category=="lcl") {
    //     keywords["color"] = "red";
    // } else if (category=="unispg") {
    //     keywords["color"] = "red";
    // }
    // keywords["hatch"] = "-";

    plt::plot(line_left_wrp, ylims_left_wrp, keywords);
    plt::plot(line_right_wrp, ylims_right_wrp, keywords);
    // plt::plot([mid, span[1]], [self.ylims_unispg["intron_max"], self.ylims_unispg["intron_min"]],
    //                     c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
    return true;
}

bool _draw_exon(map<string, float>& ylims, vector<int>& span, map<string, string>& keywords) {
    // std::map<string, string> keywords;
    // keywords["alpha"] = "0.4";
    // keywords["color"] = "grey";
    // keywords["hatch"] = "-";
    // plt::fill_between(span, span, span, keywords);
    vector<int> ylims_exon_min_wrp;
    ylims_exon_min_wrp.push_back(ylims["exon_min"]);
    vector<int> ylims_exon_max_wrp;
    ylims_exon_max_wrp.push_back(ylims["exon_max"]);


    plt::fill_between(span, ylims_exon_min_wrp, ylims_exon_max_wrp, keywords);
    // ,
                            // edgecolor=bgColor, facecolor=exonColor)
    return true;
}

void draw(vector<vector<int> >& exonIntervals, vector<vector<int> >& exonIntervals_unispg, string filename) {
    // Prepare data.
    // fprintf(stderr, "Inside draw function \n ");

    // Set the size of output image to 1200x780 pixels
    // plt::subplot(15, 1.5, 0);
    // plt::axes.set_facecolor(bgColor);
    int numExons = exonIntervals.size();
    int numExons_unispg = exonIntervals_unispg.size();

    int totalSpan = exonIntervals.back()[1] - exonIntervals[0][0];
    int totalSpan_unispg = exonIntervals_unispg.back()[1] - exonIntervals_unispg[0][0];


    float minExonLen = static_cast<float>(totalSpan)*0.005;
    float minExonLen_unispg = static_cast<float>(totalSpan_unispg)*0.005;
    fprintf(stderr, "minExonLen: %f\n", minExonLen);
    map<string, float> ylims = {{"exon_max", 4.0}, {"exon_min", 3.0}};
    map<string, float> ylims_unispg = {{"exon_max", 2.0}, {"exon_min", 1.0}};
    plt::figure_size(1200, 240);


    _set_limits(ylims, ylims_unispg);
    _transform_spans(exonIntervals, exonIntervals_unispg, minExonLen, minExonLen_unispg, totalSpan, totalSpan_unispg, numExons, numExons_unispg);

    std::map<string, string> keywords_lcl;
    keywords_lcl["color"] = "red";
    std::map<string, string> keywords_unispg;
    keywords_unispg["color"] = "green";
    std::map<string, string> keywords_black;
    keywords_black["color"] = "black";
    // if (category=="lcl") {
    //     keywords["color"] = "red";
    // } else if (category=="unispg") {
    //     keywords["color"] = "red";
    // }

    for (int i = 0; i < numExons; i++) {
        if (i > 0) {
            _draw_intron(ylims, exonIntervals[i-1][1], exonIntervals[i][0], keywords_black);
            // _draw_intron(ylims_unispg, exonIntervals_unispg[i-1][1], exonIntervals_unispg[i][0]);
        }
        // for (int j = 0; j < exonIntervals[i].size(); j++) {
        //     fprintf(stderr, "exonIntervals[i][j]: %d\n", exonIntervals[i][j]);
        // } 
        _draw_exon(ylims, exonIntervals[i], keywords_lcl);
        // _draw_exon(ylims_unispg, exonIntervals_unispg[i]);
    }
    
    for (int i = 0; i < numExons_unispg; i++) {
        if (i > 0) {
            _draw_intron(ylims_unispg, exonIntervals_unispg[i-1][1], exonIntervals_unispg[i][0], keywords_black);
        }
        // for (int j = 0; j < exonIntervals[i].size(); j++) {
        //     fprintf(stderr, "exonIntervals[i][j]: %d\n", exonIntervals[i][j]);
        // } 
        _draw_exon(ylims_unispg, exonIntervals_unispg[i], keywords_unispg);
    }
    
    // Plot line from given x and y data. Color is selected automatically.
    // plt::plot(x, y);
    // Plot a red dashed line from given x and y data.
    // plt::plot(x, w,"r--");
    // Plot a line whose name will show up as "log(x)" in the legend.
    // plt::named_plot("log(x)", x, z);
    // Set x-axis to interval [0,1000000]
    // plt::xlim(0, 1000*1000);
    // Add graph title

    keywords_black["rotation"] = "horizontal";
    keywords_black["horizontalalignment"] = "right";
    keywords_black["verticalalignment"] = "center";
    // keywords_black["rotation"] = "horizontal";

    plt::title("Local graph   vs   universal graph");
    plt::ylabel("Local graph \n\n\n\n\n\n\n Universal graph", keywords_black);
    // plt::ylabel("Universal splice graph", keywords_unispg);
    // Enable legend.
    plt::legend();
    // Save the image (file format is determined by the extension)
    plt::save(filename);
    plt::close();
}