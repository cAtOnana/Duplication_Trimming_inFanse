#include"Duplication_Trimming_inFanse.h"

const string fastqlistname = "BGI_fastqnamelist.txt";
const string fanselistname = "BGI_fansenamelist.txt";
const int batch_step = 8;
int main()
{
	ifstream fastqlistin(fastqlistname);
	if (!fastqlistin.is_open())
	{
		cout << fastqlistname << "�޷��򿪣��밴������˳�\n";
		cin.get();
		exit(EXIT_FAILURE);
	}
	vector<string> fastqlist;
	while (fastqlistin)//����fastq�ļ����б�
	{
		string temp;
		getline(fastqlistin, temp);
		if (!fastqlistin) break;
		fastqlist.emplace_back(temp);
	}

	ifstream fanselistin(fanselistname);
	if (!fanselistin.is_open())
	{
		cout << fanselistname << "��ʧ�ܣ��밴������˳�\n";
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
	for (int batch = 0; batch < fanselist.size() / 2; batch = batch + batch_step)//ÿ�δ���3�Լ�6�飬�õ�6����
	{
		//��fastq�������quality_matrix
		int upperbound = (batch + batch_step < fanselist.size() / 2) ? batch + batch_step : fanselist.size() / 2;//ѡС����Ϊ�Ͻ�����ӣ���ֹ���һ��Խ��
		for (int f = 2*batch;f<2*(upperbound);f++)
		{
			ifstream fastqin(fastqlist[f]);
			cout << "Opening " << fastqlist[f] << "...\n";
			if (!fastqin.is_open())
			{
				cout << fastqlist[f] << "�޷��򿪣��밴������˳���\n";//��Ϊfastq�����ļ���fanse�����ļ�����Ҫ�ж�Ӧ��ϵ������ĳһ�ļ��Ĵ���ʽ��������ֶ�Ӧ�������˴���
				cin.get();
				exit(EXIT_FAILURE);
			}
			vector <fastq> temp2Mfastq_reads;
			vector<quality_inform> temp_q_list;//�����洢��ǰfastq�ļ�����ȡ�����ݣ�����Ტ��quality_matrix��
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
	//��ʼ����ʹ���fanse�ļ�
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
		//�ֿ������ܣ���Ȼ�ڴ濸��ס2333
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