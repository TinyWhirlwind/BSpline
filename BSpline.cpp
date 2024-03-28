#include<iostream>
#include<cmath>
#include<vector>
#include<map>
#include<Eigen/Eigen>
#include<Eigen/SparseLU>
#include<OpenMesh\Core\IO\MeshIO.hh>
#include<OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh>
#include<OpenMesh\Core\Mesh\PolyMesh_ArrayKernelT.hh>
using namespace OpenMesh;
using namespace std;
typedef TriMesh_ArrayKernelT<> MyMesh;
///
/// de Boor-Cox递推定义
///
/// 
///节点矢量 knot
//std::map<int, double> knot;
//11个控制点 阶数k=3 节点表大小 knot.size() = 15
//Clamped列表:节点矢量 T = [0,0,0,0,1/7,2/7,3/7,4/7,5/7,6/7,1,1,1,1]
//顺序列表:节点矢量 T = [0,1/14,2/14,3/14,4/14,5/14,6/14,7/14,8/14,9/14,10/14,11/14,12/14,13/14,1]

std::vector<double> knot;
double BasisFunctionValue(int i, int k, double t)
{
	double value1, value2, value=-1.0f;//递推公式第一项和第二项的值、基函数的值
	if (k == 0)
	{
		if (t >= knot[i] && t < knot[i + 1])
			return 1;
		else
			return 0;
	}
	if (k > 0)
	{
		if (t < knot[i] || t >= knot[i + k + 1])
			return 0.0;
		float denominator = 0.0;
		float partFNum, partSNum = 0.0;//系数的值
		denominator = knot[i + k] - knot[i];//递推公式第一项分母
		if (denominator == 0.0)
		{
			partFNum = 0.0;
		}
		else
		{
			partFNum = (float)(t - knot[i]) / denominator;
		}
		denominator = knot[i + k + 1] - knot[i + 1];//递推公式第二项分母
		if (denominator == 0.0)
		{
			partSNum = 0.0;
		}
		else
		{
			partSNum = (float)(knot[i + k + 1] - t) / denominator;
		}
		value1 = partFNum * BasisFunctionValue(i, k - 1, t);
		value2 = partSNum * BasisFunctionValue(i + 1, k - 1, t);
		value = value1 + value2;
		cout << "value" << value << endl;
		return value;
	}
}
void DrawSpline(int i,int k,std::map<int,MyMesh::Point>& pointList, MyMesh& mesh)
{
	MyMesh::Point endPts = MyMesh::Point{ 0.0,0.0,0.0 };
	for (int m = k; m <= i; m++)
	{
		for (double t = knot[m]; t <= knot[m + 1]; t += 0.001)
		{
			MyMesh::Point curPts = MyMesh::Point{ 0.0,0.0,0.0 };

			for (int j = 0; j <= i; j++)
			{
				double value = BasisFunctionValue(j, k, t);
				cout << "value" << value << endl;
				curPts += (pointList[j] * value);
			}
			mesh.add_vertex(curPts);
			//endPts = curPts;
		}
		//mesh.add_vertex(endPts);
	}
}
int main()
{
	MyMesh mesh;
	int i;
	int k = 3;//3次
	//控制点集
	std::map<int, MyMesh::Point> pointList;
	pointList[0] = MyMesh::Point{ 9.29959,17.5339,25.8728 };
	pointList[1] = MyMesh::Point{ 13.2625,9.29606,22.882 };
	pointList[2] = MyMesh::Point{ 12.3762,1.28114,20.931 };
	pointList[3] = MyMesh::Point{ 6.00851,-2.94615,23.7303 };
	pointList[4] = MyMesh::Point{ -0.334078,-5.02032,24.8108 };
	pointList[5] = MyMesh::Point{ -9.86775,2.59783,24.8148 };
	pointList[6] = MyMesh::Point{ -12.7504,9.15557,24.429 };
	pointList[7] = MyMesh::Point{ -9.0734,16.5266,26.3431 };
	pointList[8] = MyMesh::Point{ -3.40426,19.974,27.9372 };
	i = pointList.size()-1;
	//节点矢量  T = [0,0,0,0,1/6,2/6,3/6,4/6,5/6,1,1,1,1] i+1个控制点 节点矢量t0----ti+k+1
	for (int n = 0; n <= k; n++)
	{
		knot.push_back(0.0f);
	}
	for (int m = 1; m < i - k + 1; m++)
	{
		knot.push_back((double)m / (i - k + 1));
	}
	for (int n = 0; n <= k; n++)
	{
		knot.push_back(1.0f);
	}

	DrawSpline(i, k, pointList, mesh);

	const char* fileout = "C:\\Users\\dell\\Desktop\\homework\\11t.obj";
	OpenMesh::IO::write_mesh(mesh, fileout);
}