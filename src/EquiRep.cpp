// Split more cycles
#include <bits/stdc++.h>
#define mst(a,b) memset((a),(b),sizeof(a))
#define psi pair <string, int>
#define pii pair <int, int>
#define mii map <int, int>
#define vi vector <int>
#define vvi vector <vi>
#define pb push_back
#define fi first
#define se second
#define N 12009 // Fitting repeat 100*50

using namespace std;

// General variables
string SG;              // Global variable of the original string
int L0;                 // L0 = SG.length();
int fa[N];              // Root of the positions([0]) in the DSU
// Edges in the dbg:
vector <pii> gp[N];
mii ee[N], ee2[N];
// Other graph related variables
pii ee_info;            // {edgecnt, edgtmax}
int in[N];              // Used in TOPOsort() and FindMMCycle() for in degree
int visit[N];           // Used for the DFS in FindCyclePrint()
int node_weight_sum[N]; // Sum of node weight on the DFS chain,
// Result variables
string LU;              // Long unit found (main part)
int LU_position;        // Starting position of LU in the repeat region
vector <pii> CG;        // Global variable for the cycle found
vector <int> ec[N];     // Equivalent class on the fly
int node_weight[N];     // Node weight of every class
double max_cycle_rate;  // Max coverage rate of the unit cycle found
int stk[N];             // Stack for DFS record
int stk_rev[N];         // Reverse of stk
int stk_i;              // stack counter
int stk2[N];            // Stack for DFS edge (for CG)
// Matrix related
int mat[2][N][N];       // Matrix values
struct Matrix_value {
    int v;
    int i, j;
};
string maj[N];          // Majority character of a class
// Alignment parameters
int match = 3, substitution = -5, gap_open = -6, gap_extend = -5;

string RegionFind(const string& seq);
psi  FindSmallUnit(const string& seq);
string FindLengthMUnit(const string& seq, int m);
int  UnitScore(const string& seq, const string& unit, int unit_position);
int  DSU_find(int x);
int  DSU_unite(int x, int y);
bool is_largest_in_window(const vvi& dp, int i, int j, int n, int w);
int  SelfLA(const string& seq);
void Matrix_refinement();
void PathDPonMat();
void MajSet();
void ConnectEdges(int n);
void Spliting(int n, int op);
void FindMMCycle(int n);
int  TOPOsort(int thres, int n);
int  FindCyclePrint(int st, int cur, int thres, int curlen, int cur_weight_sum);
int  CanTheyCombine(int i, int j);
bool cmp_mat_v(Matrix_value x, Matrix_value y);
string RegionFindChaining(const string& seq);
void Matrix_sparse_refinement();

ofstream    oF;         // Result
ifstream    iF;         // Reads

int main(int argc, const char * argv[]) {
    
        if (argc != 3) {
        cout << "Usage: ./EquiRep <input_file> <output_file_prefix>\n";
        cout << "Description:\n";
        cout << "  This program processes an input FASTA file and generates an output file.\n";
        cout << "Arguments:\n";
        cout << "  <input_file>         Path to the input FASTA file.\n";
        cout << "  <output_file_prefix> Prefix for the output FASTA file.\n";
        cout << "Example:\n";
        cout << "  ./EquiRep input.fasta output\n";
        return 1;
    }

    string fn1 = argv[1];
    ifstream iF(fn1);
    if (!iF.is_open()) {
        cerr << "Error: Unable to open input file " << fn1 << "\n";
        return 1;
    }

    string out_fn = argv[2];
    string output_fn = out_fn + ".fasta";
    ofstream oF(output_fn);
    if (!oF.is_open()) {
        cerr << "Error: Unable to create output file " << output_fn << "\n";
        return 1;
    }

    string seq;
    int unit_found_cnt = 0;
    string header;
    while (getline(iF, header)) {
        getline(iF, seq);
        SG = RegionFind(seq);
        seq = SG;
        if (seq.length() > 12000) {
            oF << ">" << header << "\n";
            oF << "toolong"<< "\n";
            continue;
        }
        if (seq.length() <= 20) {
            oF << ">" << header << "\n";
            oF << "nocycle"<< "\n";
            continue;
        }
        psi SU = FindSmallUnit(seq);
        oF << header << "\n";
        // oF << SU.fi << ' ' << SU.se << '\n';
        L0 = SG.length();
        for (int i = 0; i < L0; ++ i) {
            ec[i].clear();
            ec[i].pb(i);
            fa[i] = i;
        }
        // Self alignment to get the initial mat[][]
        mst(mat, 0);
        int f25 = SelfLA(seq);
        if (!f25) {
            oF << "nocycle" << '\n'; // lu_found
            continue;
        }
        for (int i = 0; i < L0; ++ i)
            DSU_find(i);
        // Matrix refinement
        Matrix_sparse_refinement();
        // Matrix_refinement();
        // dp for finding pathes in mat to get mat_v
        PathDPonMat();

        MajSet();
        for (int i = 0; i < L0; ++ i)
            sort(ec[i].begin(), ec[i].end());
        LU = "";
        CG.clear();
        mst(node_weight, 0);
        ConnectEdges(L0);
        FindMMCycle(L0);
        string lu0 = LU;
        int LU_score = -1;
        if (LU != "0") {
            LU_score = UnitScore(seq, LU, LU_position);
        }
        for (int i = 0; i < 5; ++ i) {
            Spliting(L0, i);
            MajSet();
            for (int i = 0; i < L0; ++ i)
                sort(ec[i].begin(), ec[i].end());
            LU = "";
            CG.clear();
            mst(node_weight, 0);
            ConnectEdges(L0);
            FindMMCycle(L0);
            if (LU != "0") {
                int tmp_score = UnitScore(seq, LU, LU_position);
                if (LU_score == -1 || tmp_score < LU_score) {
                    LU_score = tmp_score;
                    lu0 = LU;    
                }
            }
            if (SG.length() < 150)
                break;
        }
        LU = lu0;

        double final_score = 0;
        string final_unit = "";
        bool lu_found = false;
        if (LU != "0") {
            lu_found = true;
            if (LU_score >= SU.se) {
                final_unit = SU.fi;
                final_score = SU.se;
            }
            else {
                final_unit = LU;
                final_score = LU_score;
            }
        }
        else {
            final_unit = SU.fi;
            final_score = SU.se;
        }
        double eval_final_unit = final_score / SG.length();
        if (lu_found && eval_final_unit < 0.5) {

        }
        else {
            final_unit = "nocycle";
        }
        if (final_unit != "nocycle")
            ++ unit_found_cnt;
        oF << final_unit << '\n';
        // oF << "Region: " << SG.length() << "; rate0: " << eval_final_unit << '\n';
    }
    return 0;
}

/*
    Find the rpeat region
*/
string RegionFind(const string& seq) {
    int n = seq.length();

    if(n > 3100)
        return RegionFindChaining(seq);

    vvi dp(2, vi(n + 1, 0));
    vvi gap_in_seq1(2, vi(n + 1, INT_MIN));
    vvi gap_in_seq2(2, vi(n + 1, INT_MIN));
    int max_score = 0;
    int max_i = 0, max_j = 0;
    for (int i = 0; i < 2; ++ i) {
        dp[i][0] = 0;
        gap_in_seq1[i][0] = 0;
        gap_in_seq2[i][0] = 0;
    }
    for (int j = 0; j <= n; ++ j) {
        dp[0][j] = 0;
        gap_in_seq1[0][j] = 0;
        gap_in_seq2[0][j] = 0;
    }
    int ci = 1;
    for (int i = 1; i <= n; ++ i) {
        for (int j = 1; j <= n; ++ j) {
            int scoreDiag = dp[1 - ci][j - 1] + (seq[i - 1] == seq[j - 1] ? match : substitution);
            if (i >= j)
                scoreDiag = -1;   
            int scoreLeft = max({0, dp[ci][j - 1] + gap_open, gap_in_seq1[ci][j - 1] + gap_extend});
            int scoreUp = max({0, dp[1 - ci][j] + gap_open, gap_in_seq2[1 - ci][j] + gap_extend});
            dp[ci][j] = max({0, scoreDiag, scoreLeft, scoreUp});
            gap_in_seq1[ci][j] = max(dp[ci][j - 1] + gap_open, gap_in_seq1[ci][j - 1] + gap_extend);
            gap_in_seq2[ci][j] = max(dp[1 - ci][j] + gap_open, gap_in_seq2[1 - ci][j] + gap_extend);
            if (dp[ci][j] > max_score) {
                max_score = dp[ci][j];
                max_i = i;
                max_j = j;
            }
        }
        ci = 1 - ci;
    }
    string seq2 = seq.substr(0, max_j);
    reverse(seq2.begin(), seq2.end());

    n = seq2.length();
    vvi dp2(2, vector<int>(n + 1, 0));
    vvi gap_in_seq12(2, vector<int>(n + 1, INT_MIN));
    vvi gap_in_seq22(2, vector<int>(n + 1, INT_MIN));
    
    int max_score2 = 0;
    int max_i2 = 0, max_j2 = 0;
    
    for (int i = 0; i < 2; ++ i) {
        dp2[i][0] = 0;
        gap_in_seq12[i][0] = 0;
        gap_in_seq22[i][0] = 0;
    }
    for (int j = 0; j <= n; ++ j) {
        dp2[0][j] = 0;
        gap_in_seq12[0][j] = 0;
        gap_in_seq22[0][j] = 0;
    }
    ci = 1;
    for (int i = 1; i <= n; ++ i) {
        for (int j = 1; j <= n; ++ j) {
            int scoreDiag = dp2[1 - ci][j - 1] + (seq2[i - 1] == seq2[j - 1] ? match : substitution);
            if (i >= j)
                scoreDiag = -1;   
            int scoreLeft = max({0, dp2[ci][j - 1] + gap_open, gap_in_seq12[ci][j - 1] + gap_extend});
            int scoreUp = max({0, dp2[1 - ci][j] + gap_open, gap_in_seq22[1 - ci][j] + gap_extend});
            dp2[ci][j] = max({0, scoreDiag, scoreLeft, scoreUp});
            gap_in_seq12[ci][j] = max(dp2[ci][j - 1] + gap_open, gap_in_seq12[ci][j - 1] + gap_extend);
            gap_in_seq22[ci][j] = max(dp2[1 - ci][j] + gap_open, gap_in_seq22[1 - ci][j] + gap_extend);
            if (dp2[ci][j] > max_score2) {
                max_score2 = dp2[ci][j];
                max_i2 = i;
                max_j2 = j;
            }
        }
        ci = 1 - ci;
    }
    seq2 = seq2.substr(0, max_j2);
    reverse(seq2.begin(), seq2.end());
    return seq2;
}

/*
    Find small units (length 2-6)
*/
psi FindSmallUnit(const string& seq) {
    int SU_score = 0;
    string SU_min = "";
    for (int i = 2; i <= 6; ++ i) {
        string tmp_SU = FindLengthMUnit(seq, i);
        int tmp_score = UnitScore(seq, tmp_SU, 0);
        if (i == 2 || tmp_score < SU_score) {
            SU_score = tmp_score;
            SU_min = tmp_SU;
        }
    }
    return {SU_min, SU_score};
}

/*
    Find a small unit of length m
*/
string FindLengthMUnit(const string& seq, int m) {
    map <string, int> kmer_count;
    kmer_count.clear();
    for (int i = 0; i <= seq.length() - m; ++ i) {
        string tmp = seq.substr(i, m);
        string mi = tmp;
        for (int j = 1; j < m; ++ j) {
            rotate(tmp.begin(), tmp.begin() + 1, tmp.end());
            if (tmp < mi)
                mi = tmp;
        }
        ++ kmer_count[mi];
    }
    int ma = 0;
    string ma_unit;
    for (auto x : kmer_count) {
        if (x.se > ma) {
            ma = x.se;
            ma_unit = x.fi;
        }
    }
    return ma_unit;
}

/*
    Evalution of a unit.
    Edit distance between seq and concatenation of unit
*/
int UnitScore(const string& seq, const string& unit, int unit_position) {
    int n = seq.length(), n2 = unit.length(), offset1 = unit_position % n2;
    string seq2 = "";
    seq2 = unit.substr(n2 - offset1);
    for (int i = 0; i < (n - offset1) / n2; ++ i) {
        seq2 += unit;
    }
    seq2 += unit.substr(0, (n - offset1) % n2);
    int m = seq2.length();

    vvi dp(n + 1, vi(m + 1));

    // Fill dp[][] bottom up
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            // If the first string is empty, insert all characters of the second string
            if (i == 0) {
                dp[i][j] = j;
            }
            // If the second string is empty, remove all characters of the first string
            else if (j == 0) {
                dp[i][j] = i;
            }
            // If last characters are the same, ignore and recurse for the remaining strings
            else if (seq[i - 1] == seq2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            }
            // If the last character is different, consider all possibilities and find the minimum
            else {
                dp[i][j] = 1 + min({dp[i - 1][j],       // Remove
                                    dp[i][j - 1],       // Insert
                                    dp[i - 1][j - 1]}); // Replace
            }
        }
    }
    // Return the edit distance between the two strings
    return dp[m][n];
}

/*
    Disjoint set
*/
int DSU_find(int x) {
    return fa[x] == x ? x : fa[x] = DSU_find(fa[x]);
}
int DSU_unite(int x, int y) {
    int tmp1 = DSU_find(x), tmp2 = DSU_find(y), f = 0;
    if (tmp1 > tmp2)
        swap(tmp1, tmp2);
    else if (tmp1 == tmp2)
        return 0;
    if (CanTheyCombine(tmp1, tmp2) == 0)
        return 0;
    fa[tmp1] = tmp2;

    ec[tmp2].insert( ec[tmp2].end(), ec[tmp1].begin(), ec[tmp1].end() );
    ec[tmp1].clear();
    return 1;
}

/*
    Function to check if dp[i][j] is the largest in its 3x3 window
*/
bool is_largest_in_window(const vvi& dp, int i, int j, int n, int w) {
    int curr_value = dp[i][j];
    int half_w = w / 2;
    for (int di = -half_w; di <= half_w; ++ di) {
        for (int dj = -half_w; dj <= half_w; ++ dj) {
            int ni = i + di, nj = j + dj;
            if (ni > 0 && ni <= n && nj > 0 && nj <= n) { // Ensure indices are within bounds
                if (dp[ni][nj] > curr_value) {
                    return false;
                }
            }
        }
    }
    return true;
}

/*
    Self local alignment and combine position wise equiv classes
*/
int SelfLA(const string& seq) {
    int n = seq.length();
    vvi dp(n + 1, vector<int>(n + 1, 0));
    vvi gap_in_seq1(n + 1, vector<int>(n + 1, INT_MIN));
    vvi gap_in_seq2(n + 1, vector<int>(n + 1, INT_MIN));
    int max_score = 0;
    int max_i = 0, max_j = 0;
    for (int i = 0; i <= n; ++ i) {
        dp[i][0] = 0;
        gap_in_seq1[i][0] = 0;
        gap_in_seq2[i][0] = 0;
    }
    for (int j = 0; j <= n; ++ j) {
        dp[0][j] = 0;
        gap_in_seq1[0][j] = 0;
        gap_in_seq2[0][j] = 0;
    }
    for (int i = 1; i <= n; ++ i) {
        for (int j = 1; j <= n; ++ j) {
            int scoreDiag = dp[i - 1][j - 1] + (seq[i - 1] == seq[j - 1] ? match : substitution);
            if (i >= j)
                scoreDiag = -1;   
            int scoreLeft = max({0, dp[i][j - 1] + gap_open, gap_in_seq1[i][j - 1] + gap_extend});
            int scoreUp = max({0, dp[i - 1][j] + gap_open, gap_in_seq2[i - 1][j] + gap_extend});
            dp[i][j] = max({0, scoreDiag, scoreLeft, scoreUp});
            gap_in_seq1[i][j] = max(dp[i][j - 1] + gap_open, gap_in_seq1[i][j - 1] + gap_extend);
            gap_in_seq2[i][j] = max(dp[i - 1][j] + gap_open, gap_in_seq2[i - 1][j] + gap_extend);
            if (dp[i][j] > max_score) {
                max_score = dp[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
    int t25 = 25; // * n / 1000;
    // Select values in dp and sort
    vector <Matrix_value> dp_v;
    for (int i = 1; i <= n; ++ i)
        for (int j = 1; j <= n; ++ j)
            if (dp[i][j] > 0 && is_largest_in_window(dp, i, j, n, 7)) {
                Matrix_value tmp;
                tmp.i = i;
                tmp.j = j;
                tmp.v = dp[i][j];
                dp_v.pb(tmp);
            }
    sort(dp_v.begin(), dp_v.end(), cmp_mat_v);
    // Trace back and reduce path score at meeting points.
    int f25 = 0;
    for (auto x : dp_v) {
        if (x.v < t25)
            break;
        int curr_score = x.v;
        int i = x.i, j = x.j;
        vector<pii> path_positions;  // Store positions of the current path

        while (i > 0 && j > 0) {
            if (dp[i][j] == dp[i - 1][j - 1] + (seq[i - 1] == seq[j - 1] ? match : substitution)) {
                // i-1 -> j-1
                path_positions.push_back({i, j});
                // Check if meeting another path
                if (mat[0][i][j] > 0)  {
                    if(curr_score - dp[i][j] < t25) {
                        
                    }
                    else {
                        if (!f25) {
                            f25 = 1;
                        }
                        for (const auto& pos : path_positions)
                            if (mat[0][pos.fi][pos.se] == 0)
                                mat[0][pos.fi][pos.se] = curr_score - dp[i][j];
                    }
                    break;
                }
                -- i;
                -- j;
            } else if (dp[i][j] == gap_in_seq1[i][j]) {
                -- j;
            } else {
                -- i;
            }
            // Stop tracing back if reaching the start of the sequence
            if (dp[i][j] == 0) {
                if (!f25) {
                    f25 = 1;
                }
                for (const auto& pos : path_positions) 
                    if (mat[0][pos.fi][pos.se] == 0)
                        mat[0][pos.fi][pos.se] = curr_score;
                break;
            }
        }
    }
    return f25;
}

/*
    Refine the matrix mat[][]
*/
void Matrix_refinement() {
    int base_range = 100;
    int mat_max = 0;
    for (int i = 0; i < L0; ++ i)
        for (int j = 0; j < L0; ++ j)
            mat_max = max(mat_max, mat[0][i][j]);
    for (int i = 0; i < L0; ++ i)
        for (int j = 0; j < L0; ++ j)
            mat[0][i][j] = mat[0][i][j] * base_range / mat_max;

    vector <int> non_zero_entry[L0];
    for (int i = 0; i < L0; ++ i) {
        non_zero_entry[i].clear();
        for (int j = 0; j < L0; ++ j) {
            if (mat[0][i][j] > 0)
                non_zero_entry[i].pb(j);
        }
    }

    // Matrix operations
    int row_sum[N];
    for (int i = 0; i < L0; ++ i) {
        row_sum[i] = 0;
        for (int j = 0; j < L0; ++ j)
            row_sum[i] += mat[0][i][j];
    }
    int op = 0;
    for (int rep = 0; rep < 5; ++ rep) {
        for (int i = 0; i < L0; ++ i) {
            int tmp_sum = 0;
            for (int j = 0; j < L0; ++ j)
                tmp_sum += mat[op][i][j];
            for (int j = 0; j < L0; ++ j)
                if (mat[op][i][j] > 0 && tmp_sum > 0)
                    mat[op][i][j] = mat[op][i][j] * row_sum[i] / tmp_sum;
        }
        
        for (int i = 0; i < L0; ++ i)
            for (int j = 0; j < L0; ++ j)
                mat[1 - op][i][j] = mat[op][i][j];
        for (int i = 0; i < L0; ++ i) {
            for (auto x : non_zero_entry[i]) {
                for (auto y : non_zero_entry[i]) {
                    if (x == y)
                        continue;
                    if (mat[op][x][y] > 0) {
                        int tmp = min(mat[op][i][x], mat[op][i][y]);
                        tmp = min(tmp, mat[op][x][y]);
                        mat[1 - op][i][x] += tmp; 
                        mat[1 - op][i][y] += tmp;
                        mat[1 - op][x][y] += tmp;
                    }
                }
            }
        }
        op = 1 - op;
    }
    for (int i = 0; i < L0; ++ i)
        for (int j = 0; j < L0; ++ j)
            mat[1 - op][i][j] = mat[op][i][j];
    for (int i = 0; i < L0; ++ i) {
        int tmp_sum = 0;
        for (int j = 0; j < L0; ++ j)
            tmp_sum += mat[0][i][j];
        for (int j = 0; j < L0; ++ j)
            if (mat[0][i][j] > 0 && tmp_sum > 0)
                mat[0][i][j] = mat[0][i][j] * row_sum[i] / tmp_sum;
    }
    mat_max = 0;
    for (int i = 0; i < L0; ++ i)
        for (int j = 0; j < L0; ++ j)
            mat_max = max(mat_max, mat[0][i][j]);
    for (int i = 0; i < L0; ++ i)
        for (int j = 0; j < L0; ++ j)
            mat[0][i][j] = mat[0][i][j] * base_range / mat_max;
    // Matrix done
}

bool dp_mk[N][N];    // Mask of pathdp[][]
int pathdp[N][N];    // DP table for find
/*
    Use DP to find path in mat[][]
    Use threshold to get mat_v[]
*/
void PathDPonMat() {
    for (int i = 0; i < L0; ++ i)
        for (int j = 0; j < L0; ++ j)
            if (mat[0][i][j] < 10) // Threshold was 0.1
                mat[0][i][j] = 0;
    mst(pathdp, 0);
    for (int i = 0; i < L0; ++ i)
        for (int j = 0; j < L0; ++ j) {
            if (mat[0][i][j] == 0)
                continue;
            if (i == 0) {
                pathdp[i][j] = mat[0][i][j];
                continue;
            } // why this?
            pathdp[i][j] = mat[0][i][j] + pathdp[i - 1][j];
            if (j > 0) {
                pathdp[i][j] = max(pathdp[i][j], mat[0][i][j] + pathdp[i - 1][j - 1]);
                pathdp[i][j] = max(pathdp[i][j], mat[0][i][j] + pathdp[i - 1][j]);
            }
        }
    mst(dp_mk, false);
    for (int i = L0 - 1; i >= 0; -- i) {
        for (int j = L0 - 1; j >= 0; -- j) {
            if (dp_mk[i][j])
                continue;
            if (pathdp[i][j] > 100) {
                dp_mk[i][j] = true;
                int ii = i, jj = j;
                while (ii > 0 && jj >= 0) {
                    int ff = 0;
                    if (jj > 0 && !dp_mk[ii][jj - 1] && 
                        pathdp[ii][jj] == pathdp[ii][jj - 1] + mat[0][ii][jj]) {
                        -- jj;
                        dp_mk[ii][jj] = true;
                        ff = 1;
                    }
                    if (jj > 0 && ii > 0 && !dp_mk[ii - 1][jj - 1] && 
                        pathdp[ii][jj] == pathdp[ii - 1][jj - 1] + mat[0][ii][jj]) {
                        -- jj;
                        -- ii;
                        dp_mk[ii][jj] = true;
                        ff = 1;
                    }
                    if (ii > 0 && !dp_mk[ii - 1][jj] && 
                        pathdp[ii][jj] == pathdp[ii - 1][jj] + mat[0][ii][jj]) {
                        -- ii;
                        dp_mk[ii][jj] = true;
                        ff = 1;
                    }
                    if (!ff)
                        break;
                }
            }
        }
    }
    vector <Matrix_value> mat_v;
    int num = 0;
    for (int i = 0; i < L0; ++ i)
        for (int j = i + 1; j < L0; ++ j)
            if (dp_mk[i][j] == true && mat[0][i][j] > 10) { // if (pathdp[i][j] > 1) { // Threshold
                Matrix_value tmp;
                tmp.i = i;
                tmp.j = j;
                tmp.v = mat[0][i][j]; // pathdp[i][j];
                mat_v.pb(tmp);
            }
    sort(mat_v.begin(), mat_v.end(), cmp_mat_v);
    for (auto x : mat_v) {
        int i = x.i, j = x.j;
        DSU_unite(i, j);
    }
}

void MajSet() {
    for (int i = 0; i < L0; ++ i) {
        int char_cnt[30];
        mst(char_cnt, 0);
        for (auto x : ec[i]) {
            ++ char_cnt[SG[x] - 'A'];
        }
        int cnt_max = 0;
        char char_max;
        for (int j = 0; j < 26; ++ j)
            if (char_cnt[j] > cnt_max) {
                cnt_max = char_cnt[j];
                char_max = 'A' + j;
            }
        maj[i] = char_max;
    }
}

/*
    Helper functions for ConnectEdge()
    Change edges form from
        gp[] ({1, 1, 2, 2, 2, 3})
        to
        ee[] ({ {1,2}, {2,3}, {3,1} })
*/
bool GP_CMP(pii p1, pii p2) {
    return p1.fi < p2.fi;
}
bool EE_CMP(pii p1, pii p2) {
    return p1.se > p2.se;
}
void GP2EE(int n) {
    int edgecnt = 0, edgemax = 0;
    for (int i = 0; i < n; ++ i) {
        int cur = -1, curcnt = 0;
        sort(gp[i].begin(), gp[i].end(), GP_CMP);
        for (auto x : gp[i]) {
            if (cur < 0) {
                cur = x.fi;
                curcnt = x.se;
            }
            else if (cur == x.fi) {
                curcnt += x.se;
            }
            else {
                ee[i][cur] = curcnt;
                ee2[cur][i] = curcnt;
                edgemax = max(edgemax, curcnt);
                ++ edgecnt;
                cur = x.fi;
                curcnt = x.se;
            }
        }
        if (cur > 0) {
            ee[i][cur] = curcnt;
            ee2[cur][i] = curcnt;
            edgemax = max(edgemax, curcnt);
            ++ edgecnt;
        }
    }
    ee_info.fi = edgecnt;
    ee_info.se = edgemax;
}

/*
    Connect edges in dbg: (i) -> (i+1)
*/
void ConnectEdges(int n) { // n: number of nodes
    for (int i = 0; i < n; ++ i) {
        gp[i].clear();
        ee[i].clear();
        ee2[i].clear();
    }
    for (int i = 0; i < n - 1; ++ i) {
        ++ node_weight[DSU_find(i)];
        if (i < n - 2 && ec[DSU_find(i + 1)].size() == 1)
            gp[DSU_find(i)].pb({DSU_find(i + 2), 1});
        else    
            gp[DSU_find(i)].pb({DSU_find(i + 1), 1});
    }
    GP2EE(n);
}

/*
    Split nodes
*/
void NodeSplit(int cur, int med) {
    int weight_max = 0, prev = -1;
    for (auto x : ee2[cur]) {
        if (x.se == 0)
            continue;
        if (weight_max < x.se && ec[x.fi].size()) {
            weight_max = x.se;
            prev = x.fi;
        }
    }
    if (prev < 0 || ec[prev].size() <= 1) {
        return;
    }
    
    vi ec1, ec2;
    for (auto x : ec[cur]) {
        bool f_x = false;
        for (auto y : ec[prev]) {
            if (x == y + 1) {
                ec1.pb(x);
                f_x = true;
                break;
            }
            else if (y < L0 - 2 && ec[DSU_find(y + 1)].size() == 1 && x == y + 2) {
                ec1.pb(x);
                f_x = true;
                break;
            }
        }
        if (!f_x)
            ec2.pb(x);
    }
    
    int n1 = -1, n2 = -1;
    for (auto x : ec1)
        n1 = max(n1, x);
    for (auto x : ec2)
        n2 = max(n2, x);
    if (n1 < 0 || n2 < 0)
        return;
    mii ee_1, ee_2, ee2_1, ee2_2;
    for (auto x : ee2[cur]) {
        ee[x.fi].erase(cur);
        int tmpw = 0;
        for (auto y : ec[x.fi]) {
            for (auto z : ec1) {
                if (z == y + 1)
                    ++ tmpw;
                else if (y < L0 - 2 && ec[DSU_find(y + 1)].size() == 1 && z == y + 2)
                    ++ tmpw;
            }
        }
        if (tmpw > 0) {
            ee[x.fi][n1] = tmpw;
            ee2_1[x.fi] = tmpw;
        }

        tmpw = 0;
        for (auto y : ec[x.fi]) {
            for (auto z : ec2) {
                if (z == y + 1)
                    ++ tmpw;
                else if (y < L0 - 2 && ec[DSU_find(y + 1)].size() == 1 && z == y + 2)
                    ++ tmpw;
            }
        }
        if (tmpw > 0) {
            ee[x.fi][n2] = tmpw;
            ee2_2[x.fi] = tmpw;    
        }
        
    }
    for (auto x : ee[cur]) {
        ee2[x.fi].erase(cur);
        int tmpw = 0;
        for (auto z : ec[x.fi]) {
            for (auto y : ec1) {
                if (z == y + 1)
                    ++ tmpw;
                else if (y < L0 - 2 && ec[DSU_find(y + 1)].size() == 1 && z == y + 2)
                    ++ tmpw;
            }
        }
        if (tmpw > 0) {
            ee2[x.fi][n1] = tmpw;
            ee_1[x.fi] = tmpw;
        }
        tmpw = 0;
        for (auto z : ec[x.fi]) {
            for (auto y : ec2) {
                if (z == y + 1)
                    ++ tmpw;
                else if (y < L0 - 2 && ec[DSU_find(y + 1)].size() == 1 && z == y + 2)
                    ++ tmpw;
            }
        }
        if (tmpw > 0) {
            ee[x.fi][n2] = tmpw;
            ee_2[x.fi] = tmpw;
        }
    }
    ee[n1] = ee_1;
    ee[n2] = ee_2;
    ee2[n1] = ee2_1;
    ee2[n2] = ee2_2;
    ec[n1] = ec1;
    ec[n2] = ec2;
    for (auto x : ec[n1])
        fa[x] = n1;
    for (auto x : ec[n2])
        fa[x] = n2;
}
void Spliting(int n, int op) {
    vi eq_size;
    for (int i = 0; i < n; ++ i) {
        if (ec[i].size() >= 1)
            eq_size.pb(ec[i].size());
    }
    sort(eq_size.begin(), eq_size.end());
    int med = eq_size[eq_size.size() * 2 / 3 + 1], med0 = eq_size[eq_size.size() / 2];
    for (int i = 0; i < n; ++ i) {
        if (ec[i].size() > max(med0, med * 2 - op)) {
            NodeSplit(i, med);
        }
    }
}

/*
    In a graph, find the cycle which maximize the minimum edge weight in the cycle.
*/
void FindMMCycle(int n) {
    int lb = 1, rb = ee_info.se;
    while (lb <= rb) {
        int mi = (lb + rb) >> 1;
        if (TOPOsort(mi, n) > 0)
            lb = mi + 1;
        else
            rb = mi - 1;
    }
    max_cycle_rate = 0;
    while (max_cycle_rate < 0.5 && rb > 0) {
        TOPOsort(rb, n); // Necessary
        int f2 = 0;
        mst(visit, 0);
        mst(node_weight_sum, 0);
        stk_i = 0;
        for (int i = 0; i < n; ++ i)
            if (!visit[i] && in[i] > 0 ) {
                stk_i = 0;
                FindCyclePrint(i, i, rb, 2, 0);
            }
        -- rb;
    }
    if (LU == "")
        LU = "0";
}

/*
    Topological sort
        Only use the edges with weight >= thres
*/
int TOPOsort(int thres, int n) {
    vector <int> L;
    queue<int> ss;
    L.clear();
    while (!ss.empty())
        ss.pop();
    mst(in, 0);
    for (int i = 0; i < n; ++ i) {
        for (auto j : ee[i])
            if (j.se >= thres)
                ++ in[j.fi];
    }

    for (int i = 0; i < n; ++ i)
        if (in[i] == 0)
            ss.push(i);
    while (!ss.empty()) {
        int u = ss.front();
        ss.pop();
        L.pb(u);
        for (auto i : ee[u])
            if (i.se >= thres) {
                if (--in[i.fi] == 0)
                    ss.push(i.fi);
            }           
    }
    if (L.size() == n)
        return 0;
    else
        return thres;
}

/*
    Given the existence of a cycle, find it and print it
*/
int FindCyclePrint(int st, int cur, int thres, int curlen, int cur_weight_sum) {
    visit[cur] = curlen;
    stk_rev[cur] = stk_i;
    stk[stk_i ++] = cur;
    cur_weight_sum += node_weight[cur];
    node_weight_sum[cur] = cur_weight_sum;
    for (auto x : ee[cur]) {
        if (x.se >= thres) {
            if (visit[x.fi] > 1) {
                vector <pii> cg_tmp;
                cg_tmp.clear();
                double cycle_weight = cur_weight_sum - node_weight_sum[x.fi] + node_weight[x.fi];
                double cycle_rate = cycle_weight / L0;
                if (max_cycle_rate < cycle_rate) {
                    max_cycle_rate = cycle_rate;
                    string tmp = "";
                    for (int i = stk_rev[x.fi]; i < stk_i; ++ i) {
                        int tmp_edge = stk2[i];
                        if (i == stk_i - 1)
                            tmp_edge = x.se;
                        tmp += maj[stk[i]];
                        cg_tmp.pb({stk[i], tmp_edge});
                    }
                    LU = tmp;
                    LU_position = x.fi;
                    CG = cg_tmp;
                }
            }
            else if (!visit[x.fi]) {
                stk2[stk_i - 1] = x.se;
                FindCyclePrint(st, x.fi, thres, curlen + 1, cur_weight_sum);
            }
        }
    }
    visit[cur] = 1;
    -- stk_i;
    return 0;
}

/*
    Prevent to merge two kmers close to each other
*/
int CanTheyCombine(int i, int j) {
    int f = 1;
    for (auto x : ec[i])
        for (auto y : ec[j])
            if (abs(x - y) < 5) { // Hardcode
                f = 0;
                break;
            }
    return f;
}

/*
    Comparsion function for Matrix_value
*/
bool cmp_mat_v(Matrix_value x, Matrix_value y) {
    return x.v > y.v;
}

const int kmer = 8;

static inline int alphabetIndex(char c)
{
    return 3 & ((c>>2) ^ (c>>1));
}

double gamma(double l)
{
    if(l == 0.0)
        return l;
    if(l < 0.0)
        l = -l;

    return kmer * 0.01 * l + 0.5 * log2(l);
}
string RegionFindChaining(const string& seq) {

    int n = seq.length();
    int st, en;

    unordered_map<int, vector<int>> index;
    vector<double> dp;
    vector<int> g;
    vector<pair<int, int>> matches;


    int kmermax = (1<<(kmer<<1));

    int v = 0;
    for(int i = 0; i < kmer; i++)
        v = (v<<2) + alphabetIndex(seq[i]);

    index.emplace(v, std::vector<int>(1, 0));

    int m = 0;
    double max_score = -1e10;
    int max_pos = -1;

    for(int i = kmer; i < n; i++) {

        v = ((v<<2) + alphabetIndex(seq[i])) % kmermax;

        unordered_map<int, vector<int>>::const_iterator got = index.find(v);
        st = i - kmer + 1;


        if(got != index.end()) {
            for(auto pos: got->second) {

                double val = kmer;
                int last = m;
                int iter = max(0, m - 200);

                for(int j = m-1 ; j >= iter; j--)
                    if(matches[j].first < pos){

                        int alpha = min({pos - matches[j].first, kmer});
                        double beta = gamma(st - matches[j].second - (pos - matches[j].first)) + 0.1 * max(st - matches[j].second - kmer, 0);

                        if(dp[j] + alpha - beta > val)
                        {
                            val = dp[j] + alpha - beta;
                            last = j;
                        }
                }

                m++;
                dp.push_back(val);
                g.push_back(last);
                matches.push_back({pos, st});

                if(val > max_score) {
                    max_score = val;
                    max_pos = m-1;
                }
            }
            index[v].push_back(st);
        }
        else
            index.emplace(v, std::vector<int>(1, st));

    }

    if(max_pos == -1)
        return "";

    int x = max_pos;
    while(g[x] != x)    
        x = g[x];

    st = matches[x].first;
    en = matches[max_pos].second + kmer - 1;

    return seq.substr(st, en - st + 1);
}

void Matrix_sparse_refinement() {
    int base_range = 100;
    int mat_max = 0;

    vector <int> non_zero_entry[L0];
    for (int i = 0; i < L0; ++ i) {
        non_zero_entry[i].clear();
        for (int j = 0; j < L0; ++ j) {
            if (mat[0][i][j] > 0)
                non_zero_entry[i].pb(j);
        }
    }

    for (int i = 0; i < L0; ++ i)
        for (auto x : non_zero_entry[i])
            mat_max = max(mat_max, mat[0][i][x]);

    for (int i = 0; i < L0; ++ i)
        for (auto x : non_zero_entry[i])
            mat[0][i][x] = mat[0][i][x] * base_range / mat_max;



    // Matrix operations
    int row_sum[N];
    for (int i = 0; i < L0; ++ i) {
        row_sum[i] = 0;
        for (auto x : non_zero_entry[i])
            row_sum[i] += mat[0][i][x];
    }
    int op = 0;

    for (int rep = 0; rep < 5; ++ rep) {
        for (int i = 0; i < L0; ++ i) {
            int tmp_sum = 0;
            for (auto x : non_zero_entry[i])
                tmp_sum += mat[op][i][x];

            for (auto x : non_zero_entry[i])
                if (mat[op][i][x] > 0 && tmp_sum > 0)
                    mat[op][i][x] = mat[op][i][x] * row_sum[i] / tmp_sum;
        }
        
        for (int i = 0; i < L0; ++ i)
            for (auto x : non_zero_entry[i])
                mat[1 - op][i][x] = mat[op][i][x];

        for (int i = 0; i < L0; ++ i) {
            for (auto x : non_zero_entry[i]) {
                for (auto y : non_zero_entry[i]) {
                    if (x == y)
                        continue;
                    if (mat[op][x][y] > 0) {
                        int tmp = min(mat[op][i][x], mat[op][i][y]);
                        tmp = min(tmp, mat[op][x][y]);
                        mat[1 - op][i][x] += tmp; 
                        mat[1 - op][i][y] += tmp;
                        mat[1 - op][x][y] += tmp;
                    }
                }
            }
        }
        op = 1 - op;
    }

    for (int i = 0; i < L0; ++ i)
        for (auto x : non_zero_entry[i])
            mat[1 - op][i][x] = mat[op][i][x];

    for (int i = 0; i < L0; ++ i) {
        int tmp_sum = 0;
        for (auto x : non_zero_entry[i])
            tmp_sum += mat[0][i][x];
        for (auto x : non_zero_entry[i])
            if (mat[0][i][x] > 0 && tmp_sum > 0)
                mat[0][i][x] = mat[0][i][x] * row_sum[i] / tmp_sum;
    }
    mat_max = 0;
    for (int i = 0; i < L0; ++ i)
        for (auto x : non_zero_entry[i])
            mat_max = max(mat_max, mat[0][i][x]);
    for (int i = 0; i < L0; ++ i)
        for (auto x : non_zero_entry[i])
            mat[0][i][x] = mat[0][i][x] * base_range / mat_max;
    // Matrix done
}