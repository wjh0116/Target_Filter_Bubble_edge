#pragma once
#pragma once
#include <chrono>
#include <ctime>
#include <ratio>
#include <queue>
#include <ctime>
//#include "infgraph.h"
#define e exp(1)
#define c 2*(exp(1)-2)

using namespace std::chrono;

class Math {
public:
	static double log2(int n) {
		return log(n) / log(2);
	}
	static double logcnk(int n, int k) {
		double ans = 0;
		for (int i = n - k + 1; i <= n; i++)
		{
			ans += log(i);
		}
		for (int i = 1; i <= k; i++)
		{
			ans -= log(i);
		}
		return ans;
	}
};

class algo
{
private:


public:
	

	static void main(InfGraph& g, Argument& arg)
	{

		sfmt_t sfmtSeed;
		sfmt_init_gen_rand(&sfmtSeed, rand());

		double total_spread = 0;
		double total_time = 0;
		string S1_file, S2_file, Target;
		S1_file = arg.dataset + "S1" + ".txt";
		S2_file = arg.dataset + "S2" + ".txt";
		Target = arg.dataset + "Target" + ".txt";
		ifstream inFile(S1_file);//存入S1
		if (inFile.is_open())
		{
			int number;
			while (inFile >> number)
			{
				g.seedSet1.push_back(number);
				g.isS1[number]=1;
			}
			inFile.close();
		}
		else
		{
			std::cerr << "Unable to open input file!" << std::endl;
		}
		ifstream inFile1(S2_file);//存入S2
		if (inFile1.is_open())
		{
			int number;
			while (inFile1 >> number)
			{
				g.seedSet2.push_back(number);
				g.isS2[number]=1;
			}
			inFile1.close();
		}
		else
		{
			std::cerr << "Unable to open input file!" << std::endl;
		}
		g.merge_Seed_to_one(g.seedSet1,g.n);//S1=n
		g.merge_Seed_to_one(g.seedSet2,g.n+1);//S2=n+1
		//determine candidate node
		for(int i=0;i<g.n;i++)
		{
			if(!g.isS1[i]&&!g.isS2[i])
				g.canNode.push_back(i);
		}
		g.generate_random_node(arg.Qnum);//随机生成Qnum个target node 
		for(int i=0;i<g.Q.size();i++) 
			cout<<"Q: "<<g.Q[i]<<endl;
		ofstream of(arg.res, ios::app);
		g.init_hyper_graph();
		// double inf1 = g.MC_based_estimate(g.seedSet1, 10000);
		// cout << "S1's influence: " << inf1 << endl;//mento-calro的估计值
		// double inf2 = g.MC_based_estimate(g.seedSet2, 10000);
		// cout << "S2's influence: " << inf2 << endl;//mento-calro的估计值
		int ori_cnt = g.check_cnt(arg);
		cout<<"original cnt: "<<ori_cnt<<endl;
		high_resolution_clock::time_point startTime = high_resolution_clock::now();
		g.merge_Seed(g.n,g.n+1);//把两个seed node合并为一个set
		if (arg.algo == "baseline")
		{
			//determine candidate edge			
			g.get_canEdge(g.canNode);
			//select edge
			//cout<<"x";
			for(int i=0;i<arg.k;i++)
			{
				vector<pair<pair<int, int>,int>> edge_cnt;
				for(int j:g.all_seednode)
				{
					for(int x=0;x<g.canEdge[j].size();x++)
					{
						int e1 = g.canEdge[j][x];
						g.probT[j].push_back(1);   
						g.gT[j].push_back(e1);
						//cout<<"xx"<<endl;
						vector<double> p1=g.estimate_pro(g.n,g.Q,10000);
						vector<double> p2=g.estimate_pro(g.n+1,g.Q,10000);
						// for(int node:g.Q)
						// 	cout<<"node: "<<node<<" pro: "<<p1[node]<<endl;
						// for(int node:g.Q)
						// 	cout<<"node: "<<node<<" pro: "<<p2[node]<<endl;
						int cnt=0;//这条边的marginal gain
						for(int y=0;y<g.Q.size();y++)
						{
							int q=g.Q[y];
							double s1_q=p1[q];
							double s2_q=p2[q];
							//cout<<s1_q<<"  "<<s2_q<<endl;
							if(s1_q>=arg.T1&&s2_q>=arg.T2&&g.Qflag[y]==0)
								cnt++;
						}
						if(cnt!=0)
						{
							edge_cnt.push_back(make_pair(make_pair(j, e1), cnt));
							cout<<"egde: ( "<<j<<" , "<<e1<<" ) , num: "<<cnt<<endl;
						}
						g.probT[j].pop_back();
						g.gT[j].pop_back();
					}
					
				}
				if(edge_cnt.size()!=0)
				{
					std::sort(edge_cnt.begin(), edge_cnt.end(), [](const auto& a, const auto& b) 
					{
						return a.second > b.second;
					});// 按照cnt值降序排序
					for(int i=0;i<edge_cnt.size();i++)//选边不能重复
					{
						int start_node=edge_cnt[i].first.first;
						int end_node=edge_cnt[i].first.second;
						if(g.isSelect[start_node-g.n][end_node]==0)
						{
							g.probT[start_node].push_back(1);   
							g.gT[start_node].push_back(end_node);
							g.seedEdge.push_back(make_pair(start_node, end_node));
							g.isSelect[start_node-g.n][end_node]=1;
							g.check_cnt(arg);//update Qflag
							break;
						}
						else continue;
					}
				}
				else {
					cout<<"we need decrease threshold or we dont need add edge!"<<endl;
					int ran1,ran2;
					//select edge randomly
					do {
        			ran1 = sfmt_genrand_uint32(&sfmtSeed) % 2;
					ran2 = sfmt_genrand_uint32(&sfmtSeed) % g.Q.size();
    				} while (g.isSelect[ran1][g.Q[ran2]]==1);
					g.probT[ran1+g.n].push_back(1);   
					g.gT[ran1+g.n].push_back(g.Q[ran2]);
					g.seedEdge.push_back(make_pair(ran1+g.n, g.Q[ran2]));
					g.isSelect[ran1][g.Q[ran2]]=1;
				}
				
			}
		}
		//// n -> the virtual unified node for S1------ n+1 -> the virtual unified node for S2
		if (arg.algo == "RIS_based"){
			int num=0;
			vector<int>sum1; vector<int>sum2; sum1.resize(100,0); sum2.resize(100,0);
			vector<vector<int>>cov1; vector<vector<int>>cov2; 
			//cov1[q][v] denotes the cov of v for target q
			cov1.resize(100,vector<int>()); cov2.resize(100,vector<int>());
			vector<vector<bool>>edgeMark1; edgeMark1.resize(100,vector<bool>(10000,false));
			vector<vector<bool>>edgeMark2; edgeMark2.resize(100,vector<bool>(10000,false));
			vector<int>flagQ; flagQ.resize(100,0); int ori_cnt=0;//Q是否满足条件 
			//generate RR sets for each q
			for(int i=0;i<g.Q.size();i++)
			{
				vector<vector<int>> tmp_hyperG;vector<vector<int>> tmp_hyperGT;
				tmp_hyperG.resize(g.n+2, vector<int>());
				for(int j=0;j<10000;j++)
					g.generate_one_rr(j,g.Q[i],tmp_hyperG,tmp_hyperGT);
				for(int j=0;j<g.n+2;j++)//s1 and s2 need independent
				{
					if(!g.isS1[j]&&!g.isS2[j]){
						cov1[i].push_back(tmp_hyperG[j].size()); cov2[i].push_back(tmp_hyperG[j].size());
					}
					else{
						cov1[i].push_back(0);cov2[i].push_back(0);
					}
				}
				cov1[i][g.n]=0;cov2[i][g.n]=0;cov1[i][g.n+1]=0;cov2[i][g.n+1]=0;
				sum1[i]=tmp_hyperG[g.n].size(); sum2[i]=tmp_hyperG[g.n+1].size();
				//cout<<"sum1: "<<sum1[i]<<" sum2: "<<sum2[i]<<endl;
				g.QhyperG.push_back(tmp_hyperG);
				g.QhyperGT.push_back(tmp_hyperGT);
				//filter q that satisfy requirements
				double p1=(double)sum1[i]/10000.0; double p2=(double)sum2[i]/10000.0;
				if(p1>=arg.T1&&p2>=arg.T2) flagQ[i]=1,ori_cnt++,num++;
				//remove RR sets covered by S1 and S2 respectively
				g.remove_sample(g.n,edgeMark1[i],tmp_hyperG,tmp_hyperGT,cov1[i]);
				g.remove_sample(g.n+1,edgeMark2[i],tmp_hyperG,tmp_hyperGT,cov2[i]);
			}
			//select k edges S1和S2可以分开讨论
			for(int i=0;i<arg.k;i++){
				vector<pair<pair<int, int>,int>> edge_cnt;
				//first consider S1
				for(int x=0;x<g.canNode.size();x++){
					int cnt=0; int v=g.canNode[x];
					for(int j=0;j<g.Q.size();j++){
						if(flagQ[j]) continue;
						if(sum2[j]<10000*arg.T2) continue;
						if(sum1[j]+cov1[j][v]>=10000*arg.T1) cnt++;
					}
					if(cnt!=0){
						edge_cnt.push_back(make_pair(make_pair(g.n, v), cnt));
					}
				}
				//then consider S2
				for(int x=0;x<g.canNode.size();x++){
					int cnt=0; int v=g.canNode[x];
					for(int j=0;j<g.Q.size();j++){
						if(flagQ[j]) continue;
						if(sum1[j]<10000*arg.T1) continue;
						if(sum2[j]+cov2[j][v]>=10000*arg.T2) cnt++;
					}
					if(cnt!=0){
						edge_cnt.push_back(make_pair(make_pair(g.n+1, v), cnt));
					}
				}
				if(edge_cnt.size()!=0){
					std::sort(edge_cnt.begin(), edge_cnt.end(), [](const auto& a, const auto& b) 
					{
						return a.second > b.second;
					});// 按照cnt值降序排序
					int start_node=edge_cnt[0].first.first;
					int end_node=edge_cnt[0].first.second;
					//select
					if(start_node==g.n)//s1
					{
						for(int j=0;j<g.Q.size();j++)
						{
							if(flagQ[j]) continue;
							sum1[j]+=cov1[j][end_node];
							cov1[j][end_node]=0;
							g.remove_sample(end_node,edgeMark1[j],g.QhyperG[j],g.QhyperGT[j],cov1[j]);
							double p1=(double)sum1[j]/10000.0; double p2=(double)sum2[j]/10000.0;
							if(p1>=arg.T1&&p2>=arg.T2) flagQ[j]=1,num++;
							g.probT[start_node].push_back(1);   
							g.gT[start_node].push_back(end_node);
							g.seedEdge.push_back(make_pair(start_node, end_node));
						}
					}
					else//s2
					{
						for(int j=0;j<g.Q.size();j++)
						{
							if(flagQ[j]) continue;
							sum2[j]+=cov2[j][end_node];
							cov2[j][end_node]=0;
							g.remove_sample(end_node,edgeMark2[j],g.QhyperG[j],g.QhyperGT[j],cov2[j]);
							double p1=(double)sum1[j]/10000.0; double p2=(double)sum2[j]/10000.0;
							if(p1>=arg.T1&&p2>=arg.T2) flagQ[j]=1,num++;
							g.probT[start_node].push_back(1);   
							g.gT[start_node].push_back(end_node);
							g.seedEdge.push_back(make_pair(start_node, end_node));
						}
					}
				}
				else {
					cout<<"we need decrease threshold or we dont need add edge!"<<endl;
					int v;int flag;
					for(int j=0;j<g.Q.size();j++)
					{
						if(flagQ[j]) continue;
						if(sum1[j]<arg.T1*10000)
						{
							v=g.Q[j];flag=1;
							break;
						}
						if(sum2[j]<arg.T2*10000)
						{
							v=g.Q[j];flag=2;
							break;
						}
					}
					if(flag==1){
						for(int j=0;j<g.Q.size();j++)
						{
							if(flagQ[j]) continue;
							sum1[j]+=cov1[j][v];
							
							cov1[j][v]=0;
							g.remove_sample(v,edgeMark1[j],g.QhyperG[j],g.QhyperGT[j],cov1[j]);
							g.probT[g.n].push_back(1);  
							g.gT[g.n].push_back(v);
							g.seedEdge.push_back(make_pair(g.n, v));
						}
					}
					else{
						for(int j=0;j<g.Q.size();j++)
						{
							if(flagQ[j]) continue;
							sum2[j]+=cov2[j][v];
							cov2[j][v]=0;
							g.remove_sample(v,edgeMark2[j],g.QhyperG[j],g.QhyperGT[j],cov2[j]);
							g.probT[g.n+1].push_back(1);  
							g.gT[g.n+1].push_back(v);
							g.seedEdge.push_back(make_pair(g.n+1, v));
						}
					}
				}
				cout<<"satisfy num: "<<num<<endl;
				if(num==g.Q.size())
				{
					cout<<"current k is : "<<i<<endl; 
					break;
				}
			}
		}
		if (arg.algo == "degree")
		{
			g.degree_based(arg);
		}
		if (arg.algo == "to_Q")
		{
			int budget=arg.k/2;
			for(int i=0;i<g.Q.size();i++)
			{
				int v=g.Q[i];
				g.probT[g.n].push_back(1);   
				g.gT[g.n].push_back(v);
				g.probT[g.n+1].push_back(1);   
				g.gT[g.n+1].push_back(v);
				if(i==budget-1)break;
			}
		}
		if (arg.algo == "random")
		{
			int budget=arg.k/2;
			int num=0;
			vector<int>flag; flag.resize(g.n,0);
			while(num<budget){
				int ran = sfmt_genrand_uint32(&sfmtSeed) % g.canNode.size();
				int v=g.canNode[ran];
				if(flag[v]==0){
					num++;flag[v]=1;
					g.probT[g.n].push_back(1);   
					g.gT[g.n].push_back(v);
					g.probT[g.n+1].push_back(1);   
					g.gT[g.n+1].push_back(v);
				}
			}
		}
		high_resolution_clock::time_point endTime = high_resolution_clock::now();
		duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
		total_time += (double)interval.count();
		cout << "time:" << interval.count() << endl;
		for(int i=0;i<g.Qflag.size();i++) g.Qflag[i]=0;
		int after_cnt = g.check_cnt(arg);
		of << ori_cnt << "\t" <<  after_cnt  << "\t" << (double)interval.count() <<  endl;

		of.close();
	}
};
