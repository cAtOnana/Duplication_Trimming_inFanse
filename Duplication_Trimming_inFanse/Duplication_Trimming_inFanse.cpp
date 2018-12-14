#include"Duplication_Trimming_inFanse.h"
const int BGI_quality_base = 33;

ostream & operator<<(ostream & out, const fanse & f)
{
	out << f.order << tab << f.seq << tab << f.errorseq << endl;
	out << f.strand << tab << f.name << tab << f.error_num << tab << f.mappingsite << tab << f.mappingsite_num;
	return out;
}

bool operator==(const fanse_inform & fi1, const fanse_inform & fi2)
{
	if (fi1.strand == fi2.strand &&
		fi1.mappingsite == fi2.mappingsite)
		return true;
	else
		return false;
}

bool if_mean_q_less(const fanse_inform & fp1, const fanse_inform & fp2)
{
	float fp1_mean_q = fp1.q.total_qua / fp1.q.weight;
	float fp2_mean_q = fp2.q.total_qua / fp2.q.weight;
	return fp1_mean_q <= fp2_mean_q;
}

bool operator<=(const fanse_pair_inform & fp1, const fanse_pair_inform & fp2)
{
	const quality_inform& q1_1 = fp1.file1_inform.q;
	const quality_inform& q1_2 = fp1.file2_inform.q;
	const quality_inform& q2_1 = fp2.file1_inform.q;
	const quality_inform& q2_2 = fp2.file1_inform.q;
	float fp1_mean_q = (q1_1.total_qua + q1_2.total_qua) / (q1_1.weight + q1_2.weight);
	float fp2_mean_q = (q2_1.total_qua + q2_2.total_qua) / (q2_1.weight + q2_2.weight);
	return fp1_mean_q <= fp2_mean_q;
}

ifstream & operator>>(ifstream & fansein, fanse & f)
{
	fansein >> f.order;
	fansein >> f.seq >> f.errorseq;
	fansein >> f.strand >> f.name >> f.error_num;
	fansein >> f.mappingsite >> f.mappingsite_num;
	if (f.seq.find(",") != string::npos) {
		f.seq.erase(f.seq.find(","));
		f.errorseq.erase(f.errorseq.find(","));
	}
	return fansein;
}

ifstream & operator>>(ifstream & fansein, vector<fanse>& f_list)
{
	fanse temp;
	while (fansein)
	{
		fansein >> temp;
		if (!fansein) break;
		f_list.emplace_back(move(temp));
	}

	return fansein;
}

istream & operator>>(istream & fastqin, fastq & f)
{
	fastqin >> f.name1 >> f.seq >> f.name2 >> f.quality;
	return fastqin;
}

istream & operator>>(istream & fastqin, vector<fastq>& f_list)
{
	while (fastqin)
	{
		fastq temp;
		fastqin >> temp;
		if (!fastqin)
			break;
		f_list.emplace_back(temp);
	}
	return fastqin;
}

void extract_qlist(const vector<fastq>& f_list, vector<quality_inform>& q_list)
{
	for (const auto& f : f_list)
	{
		float total_logP = 0.0;
		for (int i = 0; i < f.quality.length(); i++)
		{
			total_logP += (float(f.quality[i]) - BGI_quality_base)/10;
		}
		quality_inform temp{ total_logP, f.quality.length() };
		q_list.emplace_back(temp);
	}
}

void adapt_indel(vector<fanse>& f_list)
{
	for (auto& f : f_list)
	{
		for (int i = 0; i < f.errorseq.length(); i++)
		{
			if (f.errorseq[i] == '-')
			{
				f.seq.insert(i, 1, 'I');
			}
			else if (isupper(f.errorseq[i]))
			{
				f.seq.erase(i, 1);
				f.errorseq.erase(i, 1);
				i--;
			}
		}
	}
}

void extract_inform_list(map<int, fanse_inform>& Finform_map, const vector<fanse>& f_list, const vector<quality_inform>& q_list)
{
	for (const auto& f : f_list)
	{
		fanse_inform temp;
		temp.order = f.order;
		temp.strand = f.strand;
		if (temp.strand == "R")//��Ϊ���б��ض̹�������mappingλ���Ӧ����������3�ˣ���׼ȷ��ͨ�����·�ʽת��Ϊ��5������
			temp.mappingsite = f.mappingsite + f.seq.length() - 1;
		else
			temp.mappingsite = f.mappingsite;
		temp.q = q_list[temp.order-1];
		Finform_map.emplace(temp.order, temp);
	}
}

map<int, fanse_inform> thread_task(ifstream& fansein,const vector<quality_inform>& q_list,vector<fanse>& ALL_reads)
{
	map<int, fanse_inform> Finform_map;
	while (fansein)
	{
		fanse temp;
		fansein >> temp;
		if (!fansein) break;
		ALL_reads.emplace_back(temp);
	}
	adapt_indel(ALL_reads);
	//���ǲ���GATK��׼�����������������Ƿ��indel
	extract_inform_list(Finform_map, ALL_reads, q_list);
	return move(Finform_map);
 }

void thread_output(const vector<string>& fanselist,
	int i,
	const unordered_map<int, bool>& outputlist, 
	const unordered_map<int, bool>&outputsingle,
	const vector<fanse>& ALL_reads)
{
	string outname = fanselist[i] + "NoDup";
	ofstream out(outname);

	for (const auto& f : ALL_reads)
	{
		if (outputlist.find(f.order) != outputlist.end() ||
			outputsingle.find(f.order) != outputsingle.end())
			out << f << endl;
	}
	return;
}

void main_thread(const vector<string>& fanselist,int i ,const vector<vector<quality_inform>>& quality_matrix)
{
	unordered_map<fanse_inform, int, hash<fanse_inform>, equal_to<fanse_inform>> single_reads;
	unordered_map<fanse_pair_inform, pair<int, int>, hash<fanse_pair_inform>, equal_to<fanse_pair_inform>> pair_reads;
	int f_1 = 2 * i;
	int f_2 = 2 * i + 1;
	ifstream fansein_1(fanselist[f_1]);
	cout << "Openigng " << fanselist[f_1] << "...\n";
	if (!fansein_1.is_open())
	{
		cout << fanselist[f_1] << "�޷��򿪣��밴������˳���\n";//��Ϊfastq�����ļ���fanse�����ļ�����Ҫ�ж�Ӧ��ϵ������ĳһ�ļ��Ĵ���ʽ��������ֶ�Ӧ�������˴���
		cin.get();
		exit(EXIT_FAILURE);
	}
	cout << "Open success, reading...\n";
	ifstream fansein_2(fanselist[f_2]);
	cout << "Openigng " << fanselist[f_2] << "...\n";
	if (!fansein_2.is_open())
	{
		cout << fanselist[f_2] << "�޷��򿪣��밴������˳���\n";//��Ϊfastq�����ļ���fanse�����ļ�����Ҫ�ж�Ӧ��ϵ������ĳһ�ļ��Ĵ���ʽ��������ֶ�Ӧ�������˴���
		cin.get();
		exit(EXIT_FAILURE);
	}
	cout << "Open success, reading...\n";
	vector<fanse> ALL_reads1;
	future<map<int, fanse_inform>> future_1 = async(std::launch::async, thread_task, ref(fansein_1), ref(quality_matrix[f_1]),ref(ALL_reads1));
	vector<fanse> ALL_reads2;
	future<map<int, fanse_inform>> future_2 = async(std::launch::async, thread_task, ref(fansein_2), ref(quality_matrix[f_2]),ref( ALL_reads2));
	future_1.wait();
	future_2.wait();
	map<int, fanse_inform>&& map1 = future_1.get();
	map<int, fanse_inform>&& map2 = future_2.get();
	fansein_1.close();
	fansein_2.close();
	cout << fanselist[f_1] << " & " << fanselist[f_2] << " read in completed, processing...\n";
	for (const auto& m : map1)
	{
		fanse_pair_inform fp;
		if (map2.find(m.first) != map2.end())//�ҵõ��Ļ����Ž�pair_reads��
		{
			fp.file1_inform = m.second;
			fp.file2_inform = map2[m.first];
			//��Ϣ��ȡ��ϣ�ɾ��map2�ж�Ӧ���������single read�������ɸ���
			map2.erase(m.first);
			const auto& iter_ref_pair = pair_reads.find(fp);
			if (iter_ref_pair != pair_reads.end() &&
				fp <= iter_ref_pair->first)//���ҵ���ƽ����������ԭ���Ļ�����������
				continue;
			pair_reads[fp] = make_pair(fp.file1_inform.order, fp.file2_inform.order);//�����������ŷ��ʵ��ֶ��½�/�滻��Ϣ
		}
		else if (!(single_reads.find(m.second) != single_reads.end() &&		//û���ţ����Խ�m����single read��
			if_mean_q_less(m.second, single_reads.find(m.second)->first)))	//���б���û��ƽ���������ڴ�Ԫ�ص���ʱ�������Ԫ��
			single_reads[m.second] = m.second.order;
	}
	//ѭ�������� map1���Ԫ��ȫ����������һ�أ���map2�����ǣ�map2��ʣ���Ԫ�ؾ����޷���Ե�Ԫ�أ�ͬ���ӽ�single reads
	//ע�⣬�ڴ˴�������_1��_2�������reads��������single reads�У��൱�ڽ������reads����singleend�������ݴ���
	for (const auto& m : map2)
	{
		if (!(single_reads.find(m.second) != single_reads.end() &&
			if_mean_q_less(m.second, single_reads.find(m.second)->first)))	//���б���û��ƽ���������ڴ�Ԫ�ص���ʱ�������Ԫ��
			single_reads[m.second] = m.second.order;
	}//���ˣ��������
	 //ת��pair reads��single reads
	cout << fanselist[f_1] << " & " << fanselist[f_2] << " Process ending ,outputing...\n";
	unordered_map<int, bool> outputlist1, outputlist2, outputsingle;
	for (const auto& m : pair_reads)
	{
		outputlist1.emplace(m.second.first, true);
		outputlist2.emplace(m.second.second, true);
	}
	for (const auto& r : single_reads)
	{
		outputsingle.emplace(r.second, true);
	}
	//���
	auto out_fu1 = async(std::launch::async, thread_output, ref(fanselist), f_1, ref(outputlist1), ref(outputsingle),ref(ALL_reads1));
	auto out_fu2 = async(std::launch::async, thread_output, ref(fanselist), f_2, ref(outputlist2), ref(outputsingle), ref(ALL_reads2));
	out_fu1.wait();
	out_fu2.wait();
	cout << fanselist[f_1] << " & " << fanselist[f_2] << " Output complete! Congra.\n";
	return;
}
