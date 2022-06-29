
#ifndef FASTCKY_BEAMCKYPAR_H
#define FASTCKY_BEAMCKYPAR_H

#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
//#include <stdint.h>
using namespace std;

#define MIN_CUBE_PRUNING_SIZE 20

// #define lv

#ifdef lv
  typedef int value_type;
#define VALUE_MIN std::numeric_limits<int>::lowest()
#else
  typedef double value_type;
  #define VALUE_MIN std::numeric_limits<double>::lowest()
#endif

enum Manner {
  MANNER_NONE = 0,              // 0: empty
  MANNER_H,                     // 1: hairpin candidate
  MANNER_HAIRPIN,               // 2: hairpin
  MANNER_SINGLE,                // 3: single
  MANNER_HELIX,                 // 4: helix
  MANNER_MULTI,                 // 5: multi = ..M2. [30 restriction on the left and jump on the right]
  MANNER_MULTI_eq_MULTI_plus_U, // 6: multi = multi + U
  MANNER_P_eq_MULTI,            // 7: P = (multi)
  MANNER_M2_eq_M_plus_P,        // 8: M2 = M + P
  MANNER_M_eq_M2,               // 9: M = M2
  MANNER_M_eq_M_plus_U,         // 10: M = M + U
  MANNER_M_eq_P,                // 11: M = P
  /* MANNER_C_eq_U, */
  /* MANNER_C_eq_P, */
  MANNER_C_eq_C_plus_U,     // 12: C = C + U
  MANNER_C_eq_C_plus_P,     // 13: C = C + P
};

enum BestTypes {
  TYPE_C = 0,
  TYPE_H,
  TYPE_P,
  TYPE_M,
  TYPE_Multi,
  TYPE_M2,
};

bool cmp(const std::tuple<value_type, int, int>& a, const std::tuple<value_type, int, int>& b) 
{ 
    if (std::get<0>(a) != std::get<0>(b))
        return (std::get<0>(a) < std::get<0>(b));

    else
        return  ((std::get<1>(a) - std::get<2>(a)) < (std::get<1>(b) - std::get<2>(b)));
} 

// A hash function used to hash a pair of any kind
struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const
    {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

struct ribo_state {
    int ribo_score;
    ribo_state(): ribo_score(VALUE_MIN) {};
};


struct State {
    value_type score;
    Manner manner;

    union TraceInfo {
        int split;
        struct {
            int l1;
            int l2;
        } paddings;
    };

    TraceInfo trace;

    State(): manner(MANNER_NONE), score(VALUE_MIN) {};
    State(value_type s, Manner m): score(s), manner(m) {};

    void set(value_type score_, Manner manner_) {
        score = score_; manner = manner_;
    }

    void set(value_type score_, Manner manner_, int split_) {
        score = score_; manner = manner_; trace.split = split_;
    }

    void set(value_type score_, Manner manner_, int l1_, int l2_) {
        score = score_; manner = manner_;
        trace.paddings.l1 = l1_; trace.paddings.l2 = l2_;
    }
};


class BeamCKYParser {
public:
    int beam;

    bool no_sharp_turn;
    bool is_verbose;
    bool use_constraints; // lisiz, add constraints
    bool zuker;
    int  window_size; //2 + 1 + 2 = 5 in total, 5*5 window size.
    float zuker_energy_delta;

    int jnext_org = 1000000000;


    struct DecoderResult {
        std::string structure;
        value_type score;
        unsigned long num_states;
        double time;
    };

    BeamCKYParser(int beam_size=100,
                  bool nosharpturn=true,
                  bool is_verbose=false);

    DecoderResult parse(std::string& seq, std::vector<int>* cons);

    DecoderResult parse_alifold(std::vector<std::string> & MSA, float **ribo, std::vector<std::vector<ribo_state>> & pscore, std::vector<std::vector<int>> & a2s_fast, std::vector<std::vector<int>> & s5_fast, std::vector<std::vector<int>> & s3_fast, std::vector<std::vector<int>> & SS_fast, vector<float> & smart_gap);

    void outside(std::vector<int> next_pair[]); //for zuker subopt

private:
    void get_parentheses(char* result, std::string& seq);

    std::pair<std::string, std::string> get_parentheses_outside_real_backtrace(int i, int j, State& state_beta, std::map<std::tuple<BestTypes, int, int>, std::pair<std::string, std::string> >& global_visited_outside, std::map<std::tuple<BestTypes, int, int>, std::string>& global_visited_inside, std::set<std::pair<int,int> >& window_visited);
    std::string get_parentheses_inside_real_backtrace(int i, int j, State& state, std::map<std::tuple<BestTypes, int, int>, std::string>& global_visited_inside, std::set<std::pair<int,int> >& window_visited);


    int seq_length;

    std::vector<std::unordered_map<int, State>> bestH, bestP, bestM2, bestMulti, bestM;

    std::vector<int> if_tetraloops;
    std::vector<int> if_hexaloops;
    std::vector<int> if_triloops;

    // same as bestM, but ordered
    std::vector<std::vector<std::pair<value_type, int>>> sorted_bestM;

    // hzhang: sort keys in each beam to avoid randomness
    std::vector<std::pair<int, State>> keys;

    // hzhang: sort keys in each beam to avoid randomness
    void sort_keys(std::unordered_map<int, State> &map, std::vector<std::pair<int,State>> &sorted_keys);

    void sortM(value_type threshold,
               std::unordered_map<int, State> &beamstep,
               std::vector<std::pair<value_type, int>>& sorted_stepM);

    std::vector<State> bestC;
    std::vector<int> nucs;
    std::vector<std::vector<int>> nucs_MSA;

    void prepare(unsigned len);

    void update_if_better(State &state, value_type newscore, Manner manner) {
      if (state.score < newscore)
            state.set(newscore, manner);
    };

    void update_if_better(State &state, value_type newscore, Manner manner, int split) {
        if (state.score < newscore || state.manner == MANNER_NONE)
            state.set(newscore, manner, split);
    };

    void update_if_better(State &state, value_type newscore, Manner manner, int l1, int l2) {
        if (state.score < newscore || state.manner == MANNER_NONE)
            state.set(newscore, manner, l1, l2);
    };

    value_type beam_prune(std::unordered_map<int, State>& beamstep);

    // vector to store the scores at each beam temporarily for beam pruning
    std::vector<std::pair<value_type, int>> scores;


};

#endif //FASTCKY_BEAMCKYPAR_H
