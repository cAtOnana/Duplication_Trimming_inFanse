#include"Duplication_Trimming_inFanse.h"

const string fastqlistname = "BGI_fastqnamelist.txt";
const string fanselistname = "BGI_fansenamelist.txt";
const int batch_step = 8;
int main()
{
	ifstream fastqlistin(fastqlistname);
	if (!fastqlistin.is_open())
	{
		cout << fastqlistname << "无法打开，请按任意键退出\n";
		cin.get();
		exit(EXIT_FAILURE);
	}
	vector<string> fastqlist;
	while (fastqlistin)//输入fastq文件名列表
	{
		string temp;
		getline(fastqlistin, temp);
		if (!fastqlistin) break;
		fastqlist.emplace_back(temp);
	}

	ifstream fanselistin(fanselistname);
	if (!fanselistin.is_open())
	{
		cout << fanselistname << "打开失败，请按任意键退出\n";
	}
	vector<string> fanselist;
	while (fanselistin)
	{
		string temp;
		getline(fanselistin, temp);
		if (!fanselistin) break;
		fanselist.emplace_back(temp);
	}
	vector<vector<quality_inform>> quality_matrix;
	for (int batch = 0; batch < fanselist.size() / 2; batch = batch + batch_step)//每次处理3对即6组，用到6个核
	{
		//用fastq数据填充quality_matrix
		int upperbound = (batch + batch_step < fanselist.size() / 2) ? batch + batch_step : fanselist.size() / 2;//选小的作为上界的因子，防止最后一轮越界
		for (int f = 2*batch;f<2*(upperbound);f++)
		{
			ifstream fastqin(fastqlist[f]);
			cout << "Opening " << fastqlist[f] << "...\n";
			if (!fastqin.is_open())
			{
				cout << fastqlist[f] << "无法打开，请按任意键退出。\n";//因为fastq输入文件和fanse输入文件间需要有对应关系，跳过某一文件的处理方式会打乱这种对应，故作此处理。
				cin.get();
				exit(EXIT_FAILURE);
			}
			vector <fastq> temp2Mfastq_reads;
			vector<quality_inform> temp_q_list;//用来存储当前fastq文件中提取的数据，后面会并入quality_matrix中
			cout << "Opening success, reading...\n";
			while (fastqin)
			{
				for (int i = 0; i < two_million; i++)
				{
					fastq f_temp;
					fastqin >> f_temp;
					if (!fastqin) break;
					temp2Mfastq_reads.emplace_back(f_temp);
				}
				extract_qlist(temp2Mfastq_reads, temp_q_list);
				temp2Mfastq_reads.clear();
			}
			quality_matrix.emplace_back(move(temp_q_list));
			fastqin.close();
			cout << "Read ending...\n";
		}
		cout << "fastq batch reading complete!\n";
	//开始读入和处理fanse文件
		vector<future<void>> main_future;
		for (int i = batch; i < upperbound; i++)
		{
			main_future.emplace_back(async(launch::async, main_thread, ref(fanselist), i, ref(quality_matrix)));
		}
		for (const auto& f : main_future)
		{
			f.wait();
		}

		main_future.clear();
		quality_matrix.clear();
		//分开两半跑，不然内存扛不住2333
		/*for (int i = fanselist.size() / 2 / 2; i < fanselist.size() / 2; i++)
		{
			main_future.emplace_back(async(launch::async, main_thread, ref(fanselist), i, ref(quality_matrix)));
		}
		for (const auto& f : main_future)
		{
			f.wait();
		}*/
	}
	cout << "Done\n";
	return EXIT_SUCCESS;
}
template<>
struct hash<pair<int, int>>
{
	size_t operator()(const pair<int, int>& p)const
	{
		return (p.first + p.second) / 2;
	}
};
template<>
struct equal_to<pair<int, int>>
{
	bool operator()(const pair<int, int>& p1, const pair<int, int>& p2) const
	{
		return p1.first == p2.first;
	}
};
int min()
{
	for (int i = 5; i < 4; i++)
	{
		cout << "Yes" << endl;
	}
	cin.get();
}