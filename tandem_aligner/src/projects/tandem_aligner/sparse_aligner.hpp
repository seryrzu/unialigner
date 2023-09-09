//
// Created by Andrey Bzikadze on 05/18/22.
//

// TODO: Use bridge-find instead of numbef of paths 

#pragma once

#include "cigar.hpp"
using namespace std;

namespace tandem_aligner {

class MinFreqInterval {
    int len{0};
    int fst_freq{0}, snd_freq{0};
    int fst_coord{0}, snd_coord{0};

 public:
    MinFreqInterval(int len, int fst_freq, int snd_freq,
                    int fst_coord, int snd_coord) :
        len{len}, fst_freq{fst_freq}, snd_freq{snd_freq},
        fst_coord{fst_coord}, snd_coord{snd_coord} {}

    friend bool operator<(const MinFreqInterval &lhs,
                          const MinFreqInterval &rhs);
    friend class SparseAligner;

    [[nodiscard]] double Score() const {
        return double(len)/(fst_freq*snd_freq);
    }
};

bool operator<(const MinFreqInterval &lhs, const MinFreqInterval &rhs) {
    return lhs.fst_coord < rhs.fst_coord or
        lhs.fst_coord==rhs.fst_coord and lhs.snd_coord==rhs.snd_coord;
}

class SparseAligner {
    logging::Logger &logger;
    const std::experimental::filesystem::path output_dir;
    bool no_paths; 
    bool bridges;
    static std::vector<MinFreqInterval> Cols2Vec(const MaxDisjointIntervalCollections &cols) {
        std::vector<MinFreqInterval> vec;
        for (const MaxDisjointIntervalCollection &col : cols)
            for (const auto &[_, interval] : col)
                for (const int fst_coord : interval.GetFstCoords())
                    for (const int snd_coord : interval.GetSndCoords())
                        vec.emplace_back(interval.GetLen(),
                                         col.GetFstFreq(), col.GetSndFreq(),
                                         fst_coord, snd_coord);
        return vec;
    }

    std::vector<MinFreqInterval>
    GetAlignmentVec(const std::vector<MinFreqInterval> &vec) {
        //no_paths mode
        if (no_paths or bridges){
            return GetAlignmentVecExprt(vec);
        }
        //normal mode
        if (vec.empty()) {
            return {};
        }
        std::vector<double> scores{vec.front().Score()};
        std::vector<int> backtracks{-1};

        double max_score{scores.front()};
        int argmax_score{0};
        for (int i = 1; i < vec.size(); ++i) {
            double &score = scores.emplace_back(vec[i].Score());
            int &backtrack = backtracks.emplace_back(-1);
            for (int j = i - 1; j >= 0; --j) {
                if (vec[i].fst_coord >= vec[j].fst_coord + vec[j].len and
                    vec[i].snd_coord >= vec[j].snd_coord + vec[j].len) {
                    double new_score = scores[j] + vec[i].Score();
                    if (score < new_score) {
                        score = new_score;
                        backtrack = j;
                    }
                }
            }
            if (score > max_score) {
                max_score = score;
                argmax_score = i;
            }
        }

        std::vector<MinFreqInterval> alignment_vec;
        while (argmax_score!=-1) {
            alignment_vec.emplace_back(vec[argmax_score]);
            argmax_score = backtracks[argmax_score];
        }
        std::reverse(alignment_vec.begin(), alignment_vec.end());

        logger.debug() << "Max score " << max_score <<
            "; Alignment vector " << alignment_vec.size() << "\n";

        return alignment_vec;
    }

    Cigar AlignmentVec2Cigar(std::vector<MinFreqInterval> vec,
                             const std::string &fst, const std::string &snd) {
        vec.emplace_back(0, 0, 0, fst.size(), snd.size());
        Cigar cigar;
        int left_border_fst{0}, left_border_snd{0},
            right_border_fst{0}, right_border_snd{0};
        for (int i = 0; i < vec.size() - 1; ++i) {
            const MinFreqInterval &cur = vec[i], &next = vec[i + 1];
            right_border_fst = next.fst_coord;
            right_border_snd = next.snd_coord;
            int left_fst = cur.fst_coord, left_snd = cur.snd_coord;
            while (left_fst >= left_border_fst
                and left_snd >= left_border_snd
                and fst[left_fst]==snd[left_snd]) {
                left_fst--, left_snd--;
            }
            left_fst++, left_snd++;
            int right_fst = cur.fst_coord + cur.len,
                right_snd = cur.snd_coord + cur.len;
            while (right_fst < right_border_fst
                and right_snd < right_border_snd
                and fst[right_fst]==snd[right_snd]) {
                right_fst++, right_snd++;
            }
            cigar.extend(left_fst - left_border_fst, CigarMode::D);
            cigar.extend(left_snd - left_border_snd, CigarMode::I);
            cigar.extend(right_fst - left_fst, CigarMode::M);

            left_border_fst = right_fst;
            left_border_snd = right_snd;
        }
        cigar.extend(fst.size() - left_border_fst, CigarMode::D);
        cigar.extend(snd.size() - left_border_snd, CigarMode::I);

        cigar.AssertValidity(fst, snd);
        return cigar;
    }

    [[nodiscard]] std::pair<double, double>
    GetUncovered(const Cigar &cigar,
                 const std::vector<MinFreqInterval> &intervals,
                 const std::string &fst, const std::string &snd) const {
        std::vector<int> fst_cov(fst.size() + 1), snd_cov(snd.size() + 1);
        for (const MinFreqInterval &interval : intervals) {
            fst_cov[interval.fst_coord]++;
            fst_cov[interval.fst_coord + interval.len]--;
            snd_cov[interval.snd_coord]++;
            snd_cov[interval.snd_coord + interval.len]--;
        }

        {
            int i = 0, j = 0;
            for (const CigarFragment &fragment : cigar) {
                if (fragment.mode==CigarMode::M) {
                    fst_cov[i]++, snd_cov[j]++;
                    fst_cov[i + fragment.length]--;
                    snd_cov[j + fragment.length]--;
                    i += fragment.length, j += fragment.length;
                } else if (fragment.mode==CigarMode::I) {
                    j += fragment.length;
                } else {
                    VERIFY(fragment.mode==CigarMode::D);
                    i += fragment.length;
                }
            }
        }

        for (int i = 1; i < fst_cov.size(); ++i)
            fst_cov[i] += fst_cov[i - 1];
        for (int j = 1; j < snd_cov.size(); ++j)
            snd_cov[j] += snd_cov[j - 1];

        double fst_uncovered{0}, snd_uncovered{0};
        auto extract_uncovered = [](const std::vector<int> &cov,
                                    double &uncovered) {
          for (int i = 0; i < cov.size();) {
              if (cov[i]==0) {
                  int j = i + 1;
                  while (j < cov.size() and cov[j]==0) {
                      j++;
                  }
                  uncovered += j - i;
                  i = j + 1;
              } else {
                  i++;
              }
          }
          uncovered /= cov.size();
        };
        extract_uncovered(fst_cov, fst_uncovered);
        extract_uncovered(snd_cov, snd_uncovered);
        return {fst_uncovered, snd_uncovered};
    }

    int get_no_paths(int i, const std::vector<std::vector<int>> &backtracks, std::vector<long long> &n_paths){
        
        for (int pred: backtracks[i]){
            if (pred == -1){
            n_paths[i] += (long long)1;
            }
            else if (n_paths[pred] == 0){
                n_paths[i] += (long long) get_no_paths(pred,backtracks,n_paths);
            }
            else{
            n_paths[i] += (long long) n_paths[pred];
            }
        }
        if(n_paths[i]<1 or n_paths[i]>= LONG_LONG_MAX) 
            throw std::overflow_error("Sorry the number of paths exceeded long range, can't use --no_paths with this input");
        return n_paths[i];
    }

    void dfs_subgraph(int i, const std::vector<std::vector<int>> &backtracks, std::vector<std::vector<int>> &subgraph){
        for (int pred: backtracks[i]){
            subgraph[i].emplace_back(pred);
            if (pred != -1 && subgraph[pred].size()==0){
            dfs_subgraph(pred,backtracks,subgraph);
            }
        }
    }

    std::vector<std::vector<int>> reverse_graph(const std::vector<std::vector<int>> &graph){
        std::vector<std::vector<int>> reverse_gr;

        for (int i = 0; i<graph.size();i++){
            reverse_gr.emplace_back(vector<int>{-1});
        }

        // for every (i,j) edge
        for (int i =0; i< graph.size(); ++i ){
            for (int j: graph[i]){
            // add (j,i) edge in reverse
            if (j==-1) continue;
            // logger.info()<<"adding ("<<j<<","<<i<<")\n";
            
            if (j >= reverse_gr.size()){
                
                // logger.info()<<"resizing from "<<reverse_gr.size()<<" to "<<j<<"\n";
                
                for (int fill_i= reverse_gr.size();
                fill_i<=j; fill_i++){
                reverse_gr.emplace_back(vector<int>{-1});
                }
                // logger.info()<<"resized reverse_gr to "<< reverse_gr.size()<< "\n";
            }
            if(reverse_gr[j][0]==-1){
                reverse_gr[j].clear();
            } 
            reverse_gr[j].emplace_back(i);
            // logger.info()<<"Rev edge: (" << j<<","<<i<<")\n";
            }
        }
        return reverse_gr;
    }
        
    std::vector<std::vector<int>>
    undirected_graph(const std::vector<std::vector<int>> &graph) {
        std::vector<std::vector<int>> undir_gr(graph.size()+1);
        for (int i = 0; i < graph.size(); i++) {
            for (int j : graph[i]) {
            if(j==-1)j=graph.size(); //using graph size to denote source instead of -1 for indexing 
            undir_gr[i].emplace_back(j);
            undir_gr[j].emplace_back(i);
            }
        }
        return undir_gr;
    }

    void bridgeSearch( std::vector<std::vector<int>> &undir_gr, int v,
                        int p, std::vector<int> &low, std::vector<int> &tin,
                        int &timer
                        ,std::unordered_map<int ,std::vector<int>> &bridges) {
        
        low[v] = timer;
        tin[v] = timer;
        timer++;
        for (int w : undir_gr[v]) {
            
            if (w == p)
            continue;
            else if (tin[w] == 0) {
            // unvisited = tree edge
            bridgeSearch(undir_gr, w, v, low, tin, timer
                        , bridges
            );
            low[v] = min(low[v], low[w]);
            if (low[w] > tin[v]) {
                bridges[(v!=undir_gr.size()-1)?v:-1].emplace_back((w!=undir_gr.size()-1)?w:-1);
            }
            } else {
            // visited = back edge
            low[v] = min(low[v], tin[w]);
            }
        }
    }
    
    std::unordered_map<int, std::vector<int>> getBridges(const std::vector<std::vector<int>> backtracks) {
        std::vector<std::vector<int>>  undir_gr = undirected_graph(backtracks);
        std::vector<int> tin(undir_gr.size());
        std::vector<int> low(undir_gr.size());
        int timer = 1;
        std::unordered_map<int, std::vector<int>> bridges;
        for (int i=0;i< undir_gr.size();i++) {
            if (tin[i] == 0) {
            bridgeSearch(undir_gr, i, -2, low, tin, timer,bridges);
            }
        }
        return bridges;
    }

    /**
     * Record number of paths for kmers in CSV
    */
    void export_no_paths(const std::vector<MinFreqInterval> &vec,std::vector<long long> n_paths_fwd,std::vector<long long> n_paths_rev, string out_path){
        ofstream no_paths_file;
        no_paths_file.open(out_path, ios::app);
        no_paths_file<<"start_1"<<","<<"start_2"<<","<<"len"<<","<<"n_opt_paths_fwd"<<","<<"n_opt_paths_rev"<<","<<"n_count"<<","<<"m_count"<<"\n";

        for (int i = 0; i< vec.size(); i++){
                no_paths_file<<vec[i].fst_coord<<","<<vec[i].snd_coord<<","<<vec[i].len<<","<<n_paths_fwd[i]<<","<<n_paths_rev[i]<<","<<vec[i].fst_freq<<","<<vec[i].snd_freq<<"\n";
        }
        no_paths_file.close();
    }

    void export_bridges(std::unordered_map<int, std::vector<int>> &bridges,const std::vector<MinFreqInterval> &vec, string out_path) {
        // write the coordinate of the bridge edge
        ofstream bridges_file;
        bridges_file.open(out_path, ios::app);

        bridges_file << "#Bridges in alignment step\n";
        bridges_file<<"start_1"<<","<<"start_2"<<","<<"end_1"<<","<<"end_2"<<"\n";
        for (auto elem : bridges) {
            for (int j =0;j< elem.second.size();++j) {
                if(elem.first!=-1 && elem.second[j]!=-1){
                    // left segment index
                    int first = (vec[elem.first].fst_coord < vec[elem.second[j]].fst_coord)?elem.first:elem.second[j];
                    // right segment index
                    int second = (vec[elem.first].fst_coord > vec[elem.second[j]].fst_coord)?elem.first:elem.second[j];
                    int st1 = vec[first].fst_coord +vec[first].len;
                    int st2 = vec[first].snd_coord +vec[first].len;
                    int end1 = vec[second].fst_coord;
                    int end2 = vec[second].snd_coord;
                    bridges_file<<st1<<","<<st2<<","<<end1<<","<<end2<<"\n";
                }
                else{
                    int second = max(elem.first,elem.second[j]); // the segment that is not -1
                    int st1 = 0;
                    int st2 = 0;
                    int end1 = vec[second].fst_coord;
                    int end2 = vec[second].snd_coord;
                    bridges_file<<st1<<","<<st2<<","<<end1<<","<<end2<<"\n";
                }
            }
        }
    }

 public:

    std::vector<MinFreqInterval>
    GetAlignmentVecExprt(const std::vector<MinFreqInterval> &vec) {
                if (vec.empty()) {
                    return {};
                }
                std::vector<double> scores{vec.front().Score()};
                std::vector<int> backtracks{-1};
                std::vector<std::vector<int>> backtrack_all{{-1}};

                double max_score{scores.front()};
                std::vector<int> argmax_score_list{0};
                for (int i = 1; i < vec.size(); ++i) {
                    double &score = scores.emplace_back(vec[i].Score());
                    int &backtrack = backtracks.emplace_back(-1);
                    std::vector<int> bt_all_row{-1};

                    for (int j = i - 1; j >= 0; --j) {

                        if (vec[i].fst_coord >= vec[j].fst_coord + vec[j].len and
                            vec[i].snd_coord >= vec[j].snd_coord + vec[j].len) {
                            double new_score = scores[j] + vec[i].Score();

                            if (score < new_score) {
                                score = new_score;
                                backtrack = j;
                                // reset bt_all_row with new best edge
                                bt_all_row.clear();
                                bt_all_row.emplace_back(j);
                            }
                            else if (score==new_score){
                                // add tying edge to bt_all_row
                                bt_all_row.emplace_back(j);
                            }
                        }
                    }

                    backtrack_all.emplace_back(bt_all_row);

                    if (score > max_score) {
                        max_score = score;
                        argmax_score_list.clear();
                        argmax_score_list.emplace_back(i);
                    }
                    else if(score == max_score){
                        argmax_score_list.emplace_back(i);
                    }
                }
                int argmax_score = argmax_score_list[0];

                //TOREMOVE
                logger.info()<< "Found " << argmax_score_list.size()<<" optimum path ends\n";

                std::vector<MinFreqInterval> alignment_vec;
                while (argmax_score!=-1) {
                    alignment_vec.emplace_back(vec[argmax_score]);
                    argmax_score = backtracks[argmax_score];
                }
                std::reverse(alignment_vec.begin(), alignment_vec.end());
                
                logger.info() << "Max score " << max_score <<
                    "; Alignment vector " << alignment_vec.size() << "\n";
                
                if(bridges){
                    // find optimum subgraph
                    std::vector<std::vector<int>> backtrack_all_optimum(backtrack_all.size());
                    for (int sink: argmax_score_list){
                        dfs_subgraph(sink,backtrack_all, backtrack_all_optimum);
                    }

                    //bridge-search on optimum subgraph
                    std::unordered_map<int, std::vector<int>> bridges = getBridges(backtrack_all_optimum);
                    export_bridges(bridges,vec,output_dir/"bridges.txt");
                }

                if (no_paths){
                std::vector<long long> n_paths_fwd= std::vector<long long>(vec.size());
                for(int sink: argmax_score_list){
                    
                    get_no_paths(sink, backtrack_all, n_paths_fwd);

                    // //TOREMOVE
                    logger.info()<<n_paths_fwd[sink]<<" Optimum paths from "<<sink<<"\n"; 

                }

                // find startpoints of optimum paths
                std::vector<int> optimum_starts;
                for(int i =0; i< vec.size(); i++){
                    if (n_paths_fwd[i]!=0 && backtrack_all[i][0]==-1){

                        //TOREMOVE
                        logger.info()<<"optimum start at " <<i<<"\n";
                        optimum_starts.emplace_back(i);
                    }
                }

                std::vector<std::vector<int>> backtrack_all_optimum;
                for (int i =0; i<vec.size(); i++){
                    backtrack_all_optimum.emplace_back(std::vector<int>());
                    if(n_paths_fwd[i]>0){
                        backtrack_all_optimum[i] = backtrack_all[i];
                    }
                }


                std::vector<long long> n_paths_rev= std::vector<long long>(vec.size());
                std::vector<std::vector<int>> rev_backtrack_all_optimum = reverse_graph(backtrack_all_optimum);
                for(int start: optimum_starts){
                    get_no_paths(start, rev_backtrack_all_optimum, n_paths_rev);

                    //TOREMOVE
                    logger.info()<<"From "<<start<<" found "<<n_paths_rev[start]<<" paths\n";
                }

                //TODO: Change this
                export_no_paths(vec, n_paths_fwd, n_paths_rev, output_dir/"no_paths.csv");
                }
                return alignment_vec;
    }

    explicit SparseAligner(logging::Logger &logger, string output_dir, bool no_paths, bool bridges) : logger{logger}, output_dir{output_dir}, no_paths{no_paths}, bridges{bridges} {}

    Cigar Align(const MaxDisjointIntervalCollections &cols,
                const std::string &fst, const std::string &snd) {
        std::vector<MinFreqInterval> vec = Cols2Vec(cols);
        std::sort(vec.begin(), vec.end());
        logger.trace() << "Running quadratic on " << vec.size() << "...\n";
        const std::vector<MinFreqInterval> alignment_vec = GetAlignmentVec(vec);

        Cigar cigar = AlignmentVec2Cigar(alignment_vec, fst, snd);

        // const auto[fst_uncovered, snd_uncovered] =
        // GetUncovered(cigar, vec, fst, snd);
        // logger.info() << "Uncovered fst = " << fst_uncovered <<
        //               "; snd = " << snd_uncovered << "\n";

        // double mutability = cigar.GetMutability();
        // logger.info() << "Mutability = " << mutability << "\n";
        return cigar;
    }
};

} // namespace tandem_aligner