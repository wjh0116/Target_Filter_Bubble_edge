#include "iheap.h"
#include <queue>	//priority_queue
#include <utility>  // pair
#include <numeric>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <random>
class InfGraph: public Graph
{
private:
	//vector<bool> activated;
    vector<bool> visit;
    vector<int> visit_mark;
	
public:
    vector<vector<int>> hyperG;
    vector<vector<int>> hyperGT;
	vector<vector<int>> SampleG;
	vector<vector<int>> SampleGT;
	vector<vector<int>> hyperG_2;
	vector<vector<int>> hyperGT_2;
	vector<vector<vector<int>>>QhyperG;
	vector<vector<vector<int>>>QhyperGT;
	vector<int> canNode;//candidate node
	sfmt_t sfmtSeed;
	vector<int>Qflag;//判断q是否满足条件，如果是，设置为1，else 0
	vector<int> seedSet1;//S1
	vector<int> seedSet2;//S2
	vector<int> isS1;//这个节点是否为S1
	vector<int> isS2;//这个节点是否为S2
	vector<int> Q;//目标节点
	vector<int>all_seednode;
	vector<pair<int, int>> seedEdge;
	vector<vector<int>> isSelect;
	//vector<pair<pair<int, int>,int>> edge_cnt;//每条边对应的cnt
	//vector<int> edge_cnt;//每个candidate edge的marginal gain
    InfGraph(string folder, string graph_file): Graph(folder, graph_file)
    {
        srand(time(NULL));
        sfmt_init_gen_rand(&sfmtSeed , rand());		
        visit = vector<bool> (n+1);
		isSelect=vector<vector<int>>(2, std::vector<int>(n, 0));
		//QhyperG=vector<vector<vector<int>>>(100,vector<vector<int>>(n+2, vector<int>()));//
		//QhyperGT=vector<vector<vector<int>>>(100,vector<vector<int>>());//
        visit_mark = vector<int> (n+1);
		isS1 = vector<int> (n, 0);
		isS2 = vector<int> (n, 0);
		//edge_cnt = vector<int> (2*n, 0);//假设S1和S2各一个
		hyperG.resize(n+1, vector<int>());
		Qflag.resize(500,0);
		hyperG_2.resize(n + 1, vector<int>());
		SampleG.resize(n+1, vector<int>());
    }

    void init_hyper_graph(){
		for (auto& hyper : hyperG)hyper.clear();
		for (auto& hyperT : hyperGT)vector<int>().swap(hyperT);
		hyperGT.clear();
		for (auto& hyper : SampleG)hyper.clear();
		for (auto& hyperT : SampleGT)vector<int>().swap(hyperT);
		SampleGT.clear();

		for (auto& hyper : hyperG_2)hyper.clear();
		for (auto& hyperT : hyperGT_2)vector<int>().swap(hyperT);
		hyperGT_2.clear();
    }

	char* map_file(const char* fname, size_t& length)
	{
		int fd = open(fname, O_RDONLY);
		if (fd == -1)
			handle_error("open");

		// obtain file size
		struct stat sb;
		if (fstat(fd, &sb) == -1)
			handle_error("fstat");

		length = sb.st_size;

		char* addr = static_cast<char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
		if (addr == MAP_FAILED)
			handle_error("mmap");

		// TODO close fd at some point in time, call munmap(...)
		close(fd);
		return addr;
	}
	
	double MC_based_estimate(vector<int> rumorNode, int count)
	{
		vector<int> vis;
		vis.resize(n);
		double sum = 0;
		for (int tag = 1; tag <= count; tag++)
		{
			//vector<vector<int>> eSample;
			//eSample.resize(n);
			queue<int> q;
			for (int x : rumorNode)
				q.push(x), vis[x] = -tag, sum += 1.0;
			while (!q.empty())
			{
				int x = q.front();
				q.pop();
				for (int i = 0; i < gT[x].size(); i++)
				{
					if (probT[x][i] == 0)
						continue;
					if (probT[x][i] >= sfmt_genrand_real1(&sfmtSeed))
					{
						//eSample[x].push_back(gT[x][i]);
						if (vis[gT[x][i]] != -tag)
							q.push(gT[x][i]), vis[gT[x][i]] = -tag, sum += 1.0;
					}
				}
			}
		}
		return sum / ((double)count);
	}
	static inline double logcnk(const size_t n, size_t k)
	{
		k = k < n - k ? k : n - k;
		double res = 0;
		for (auto i = 1; i <= k; i++) res += log(double(n - k + i) / i);
		return res;
	}
	static inline double pow2(const double t)
	{
		return t * t;
	}
	void generate_random_node(int num)
	{

		// 生成1到n的序列
		std::vector<int> numbers(n);
		std::iota(numbers.begin(), numbers.end(), 1);
		
    	
		// 打乱序列
		std::mt19937 fixedGen(123); // 使用特定的种子值初始化 mt19937
    	std::shuffle(canNode.begin(), canNode.end(), fixedGen);
		//std::shuffle(numbers.begin(), numbers.end(), gen);
		int i=0;
		while (i < num) 
		{
			Q.push_back(canNode[i]);
			i++;
    	}
	}
	void merge_Seed(int S1,int S2)//将S1和S2中的节点合并为一个新的seed set
	{
		all_seednode.push_back(S1);
		all_seednode.push_back(S2);
	}
	void merge_Seed_to_one(vector<int> Seed,int node)//将S1或者S2中的seed变成一个节点
	{
		for(int S:Seed)
		{
			gT[node].push_back(S);
			gT_reverse[S].push_back(node);
			probT[node].push_back(1);
			probT_reverse[S].push_back(1);
		}
	}
	int check_cnt(Argument& arg)
	{
		vector<double> x1=estimate_pro(n,Q,10000);
		vector<double> x2=estimate_pro(n+1,Q,10000);
		int cnt=0;
		for(int y=0;y<Q.size();y++)
		{
			int q=Q[y];
			double s1_q=x1[q];
			double s2_q=x2[q];
			cout<<"s1 to q: "<<s1_q<<", s2 to q: "<<s2_q<<endl;
			if(s1_q>=arg.T1&&s2_q>=arg.T2&&!Qflag[y])
				cnt++, Qflag[y]=1;
		}
		return cnt;
	}
	void get_canEdge(vector<int> canNode)
	{
		for(int i:all_seednode)
		{
			for(int node:canNode)
			{
				canEdge[i].push_back(node);
			}
		}
	}
	///  this function is to generate one rr set
	void generate_one_rr(int hyperId, int uStart,vector<vector<int>> &tmp_hyperG,vector<vector<int>> &tmp_hyperGT)
	{
		vector<int> visit_mark; vector<bool> visit; visit = vector<bool>(n); visit_mark = vector<int>(n);
		unsigned int n_visit_mark = 0, curIdx = 0; 
		visit_mark[n_visit_mark++] = uStart; 
		visit[uStart] = true; 
		tmp_hyperG[uStart].push_back(hyperId); 
		while (curIdx < n_visit_mark) {
				int i = visit_mark[curIdx++]; 
				for (int j = 0; j < (int)gT_reverse[i].size(); j++) { 
					int v = gT_reverse[i][j]; 
					if (visit[v])continue; 
					double randDouble = sfmt_genrand_real1(&sfmtSeed); 
					if (randDouble > probT_reverse[i][j])continue; 
					visit[v] = true; 
					visit_mark[n_visit_mark++] = v; 
					tmp_hyperG[v].push_back(hyperId); 
				}
			}
		tmp_hyperGT.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark));
	}
	void remove_sample(int seed, vector<bool> &isDelete,vector<vector<int>> tmp_hyperG,vector<vector<int>> tmp_hyperGT,vector<int>&coverage)
	{
		for (auto edgeIdx : tmp_hyperG[seed]) {
				if (tmp_hyperGT[edgeIdx].size() == 0 || isDelete[edgeIdx]) continue; 
				isDelete[edgeIdx] = true; 
				for (auto nodeIdx : tmp_hyperGT[edgeIdx]) { 
					if (coverage[nodeIdx] == 0) continue;
					coverage[nodeIdx]--;
				}
			}
	}
	struct Node {
		int a;
		int b;

		Node(int a, int b) : a(a), b(b) {}
	};

	struct Compare {
		bool operator()(const Node& n1, const Node& n2) {
			return n1.b < n2.b; // ����b��ֵ���н�������
		}
	};
	void degree_based(const Argument& arg)
	{
		int budget=arg.k/2;int num=0;
		priority_queue<Node, std::vector<Node>, Compare> pq;
		for (int i = 0; i < canNode.size(); i++) {
			int v=canNode[i];
			pq.push(Node(v, deg[v]));
		}
		while(num<budget)
		{
			int v = pq.top().a;
			cout<<"deg: "<<deg[v]<<endl;
			pq.pop();
			probT[n].push_back(1);   
			gT[n].push_back(v);
			probT[n+1].push_back(1);   
			gT[n+1].push_back(v);
			num++;
		}
		
	}
	vector<double> estimate_pro(int seedset, vector<int> Target_nodes, double count)
	{
		vector<double> pro;
		pro.resize(n+2,0);//如果不是target node，都是0
		vector<int> vis;
		vis.resize(n+2);
		vector<int> sum;//如果不是target node，都是0
		sum.resize(n+2,0);
		vector<int> isTarget;
		isTarget.resize(n+2,0);
		for(int node:Target_nodes)
		{
			isTarget[node]=1;
		}
		
		for (int tag = 1; tag <= count; tag++)
		{
			queue<int> q;
			q.push(seedset), vis[seedset] = -tag;
			while (!q.empty())
			{
				int x = q.front();
				q.pop();
				for (int i = 0; i < gT[x].size(); i++)
				{
					if (probT[x][i] == 0)
						continue;
					if (probT[x][i] >= sfmt_genrand_real1(&sfmtSeed))
					{
						if (vis[gT[x][i]] != -tag)
						{
							q.push(gT[x][i]), vis[gT[x][i]] = -tag;
							if(isTarget[gT[x][i]])
								sum[gT[x][i]]++;
						}
					}
				}
			}
		}
		for(int node:Target_nodes)
		{
			//cout<<"node: "<<node<<" pro: "<<sum[node]<<endl;
			pro[node]=sum[node]/count*1.0;
		}
		return pro;
	}
};