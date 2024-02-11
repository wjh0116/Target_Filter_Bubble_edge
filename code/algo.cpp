#define HEAD_INFO

#include "sfmt/SFMT.h"
#include "head.h"

class Argument {
public:
    unsigned int k;
    unsigned int Qnum;
    string dataset;
    string res;
    string algo;
    double T1;
    double T2;
};

#include "graph.h"
#include "infgraph.h"
#include "algo.h"

void OutputSeedSetToFile(vector<int> seed_set, const Argument& arg)
{
    string seedfile = "results/res_" + arg.dataset;
    ofstream of(seedfile, ios::app);
    //of.open(seedfile);
    for (int seed : seed_set)
    {
        of << seed << endl;
    }
    of << endl;
    of.close();
}

void run_with_parameter(InfGraph& g,  Argument& arg)
{
    algo::main(g, arg);

    //OutputSeedSetToFile(g.seedSet, arg);
}
void Run(int argn, char** argv)
{
    Argument arg;

    arg.k = 0;

    for (int i = 0; i < argn; i++)
    {
        if (argv[i] == string("-k"))
            arg.k = atoi(argv[i + 1]);
        if (argv[i] == string("-Qnum"))
            arg.Qnum = atoi(argv[i + 1]);
        if (argv[i] == string("-dataset"))
            arg.dataset = argv[i + 1];
        if (argv[i] == string("-algo"))
            arg.algo = argv[i + 1];
        if (argv[i] == string("-T1"))
            arg.T1 = stod(argv[i + 1]);
        if (argv[i] == string("-T2"))
            arg.T2 = stod(argv[i + 1]);
    }
    //cout<<"1"<<arg.dataset<<endl;
    ASSERT(arg.dataset != "");
    string temp_name = arg.dataset.substr(arg.dataset.find_last_of("/") + 1);
    arg.res = "results/res_" + temp_name + "_K=" + to_string(arg.k) +"_|Q|=" + to_string(arg.Qnum) +
    "_T=" + to_string(arg.T1) + "_algo=" + arg.algo;
    arg.dataset = arg.dataset + "/";
    string graph_file;
    graph_file = arg.dataset + "graph_ic.inf";
    InfGraph g(arg.dataset, graph_file);
    run_with_parameter(g, arg);
}


int main(int argn, char** argv)
{
    __head_version = "v1";
    OutputInfo info(argn, argv);

    Run(argn, argv);
}

