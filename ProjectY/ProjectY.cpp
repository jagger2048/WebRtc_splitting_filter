// ProjectY.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <string>
#include <map>
#include<typeinfo>
using namespace std;
size_t up() {
	static int upCouter = 0;



	upCouter = 666;
	return upCouter++;
}
struct cTest {
	int a = 10;
	int b = 10;
	int c[3] = {1,2,3};
	int d[3] = { 4,5,6 };
	float * pp;
};
void test(const float*p){
	cout << sizeof(*p) << " == " << sizeof(p)<<endl;
	cout << *p << " == " << sizeof(p)<<endl;
}
int main()
{
	// 测试 vector 的构造函数

	float a[10] = {1,2,3,4,5,6,7,8,9,0};
	const float * p = a;
	vector<float> b(p, p+5);
	for (size_t i = 0; i < b.size(); i++)
	{
		cout << b.at(i)<<endl;
	}
	cout << sizeof(a) << endl;
	cout << sizeof(*p) << endl;
	cout << "next\n";
	cTest *h = (cTest*)malloc(sizeof(cTest));
	const float* pp = (float*)malloc(sizeof(float) * 6);
	cout << sizeof(pp) << " == " << sizeof(*pp) << endl;

	cout << "test2\n";
	test(a);
	free(h);
	return 0;
}