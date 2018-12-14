#pragma once
#include<string>
#include<iostream>
#include<fstream>
#include<unordered_map>
#include<map>
#include<vector>
#include<utility>
#include<map>
#include<future>
using namespace std;
const char tab = '\t';
const int two_million = 2000000;
struct fanse {
	int order;
	std::string seq;
	std::string errorseq;
	string strand;
	string name;
	int error_num;
	int mappingsite;
	int mappingsite_num;
	friend ostream& operator<<(ostream& out, const fanse& f);
};
struct quality_inform
{
	float total_qua;
	int weight; 
	
};
struct fanse_inform
{
	int order;
	string strand;
	int mappingsite;
	quality_inform q;
	friend bool operator==(const fanse_inform& fi1,const fanse_inform& fi2);
	friend bool if_mean_q_less(const fanse_inform& fp1, const fanse_inform& fp2);
};
struct fanse_pair_inform
{
	fanse_inform file1_inform;
	fanse_inform file2_inform;
	friend bool operator<=(const fanse_pair_inform& fp1, const fanse_pair_inform& fp2);//用于比较加权平均的q
};
struct fastq
{
	string name1;
	string seq;
	string name2;
	string quality;
};

ifstream& operator>>(ifstream& fansein, fanse& f);
ifstream& operator>>(ifstream& fansein, vector<fanse>& f_list);
istream& operator>>(istream& fastqin, fastq&f);
istream& operator>>(istream& fastqin, vector<fastq>& f_list);
void extract_qlist(const vector<fastq>& f_list, vector<quality_inform>& q_list);
void adapt_indel(vector<fanse>& f_list);
void extract_inform_list(map<int, fanse_inform>& Finform_map, const vector<fanse>& f_list, const vector<quality_inform>& q_list);
map<int, fanse_inform> thread_task(ifstream& fansein, const vector<quality_inform>& q_list, vector<fanse>& ALL_reads);
void thread_output(const vector<string>& fanselist,	int i,const unordered_map<int, bool>& outputlist,const unordered_map<int, bool>&outputsingle,const vector<fanse>& ALL_reads);
void main_thread(const vector<string>& fanselist, int i, const vector<vector<quality_inform>>& quality_matrix);


template<>
struct hash<fanse_inform>
{
	size_t operator()(const fanse_inform & f) const
	{
		return f.mappingsite;
	}
};
template<> 
struct equal_to<fanse_inform>
{
	bool operator()(const fanse_inform& L_key, const fanse_inform& R_key)const
	{
		if (L_key.strand == R_key.strand &&
			L_key.mappingsite == R_key.mappingsite)
			return true;
		else
			return false;
	}
};
template<>
struct hash<fanse_pair_inform>
{
	size_t operator()(const fanse_pair_inform & f) const
	{
		return (f.file1_inform.mappingsite+f.file2_inform.mappingsite)/2;
	}
};
template<>
struct equal_to<fanse_pair_inform>
{
	bool operator()(const fanse_pair_inform& L_key, const fanse_pair_inform& R_key)const
	{
		if (L_key.file1_inform == R_key.file1_inform &&
			L_key.file2_inform == R_key.file2_inform)
			return true;
		else if (L_key.file1_inform == R_key.file2_inform &&
			L_key.file2_inform == R_key.file1_inform)
			return true;
		else
			return false;
	}
};

inline map<int,fanse_inform>&& add(ifstream& a, const vector<quality_inform>& b) { return map<int, fanse_inform>(); }