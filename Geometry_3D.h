//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Geometry_3D.h
//OBJECTIVE:	The definitions of point, line and plane in 3D
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#ifndef GEOMETRY_3D_H
#define GEOMETRY_3D_H

#include<cmath>
#include<stdlib.h>
#include<vector>
#include "Hns.h"
using namespace hns;

#include "MathMatrix.h"
//---------------------------------------------------------------------------
//������ά�ռ����
class Point_3D
{
	public:

		//���ݳ�Ա
		double x, y, z;
		int flag;

		//���캯��
		Point_3D(){};
		Point_3D( double px, double py, double pz );

		//��Ա����
        Point_3D operator+( Point_3D &pt );
        Point_3D operator+( const Point_3D &pt );
        Point_3D operator+( double d );
        Point_3D operator-( double d );
        Point_3D operator-( Point_3D &pt );
        Point_3D operator-( const Point_3D &pt );
        Point_3D operator*( double d );
        Point_3D operator/( double d );
		bool operator==( Point_3D &pt );
		bool operator!=( Point_3D &pt );
        double distance_to(const Point_3D &pt)const;
		double distance_to(const double &px, const double &py,  const double &pz)const;
};

//---------------------------------------------------------------------------
//������ά��(��)��
class Line_3D
{
	public:

		//���ݳ�Ա
		 Point_3D point[2];	//��ʾ�߶ε������˵�����
		 double len;					//��ʾ�߶εĳ���
		 bool virtual_line;		//��ʾ�߶��Ƿ��϶�Ϊ���߶�(false:���߶Σ�����Ϊһ�㣩; true:�����߶�)


		//���캯��
		Line_3D(Point_3D p0, Point_3D p1);
		
		//��Ա����
		double length();	//�߶εĳ���
        double distance_point_to_line(const Point_3D *point_temp)const;    //�㵽ֱ�ߵľ���
        double distance_point_to_line(const Point_3D &point_temp)const;    //�㵽ֱ�ߵľ���
        double distance_point_to_line(const double dx, const double dy, const double dz)const;  //�㵽ֱ�ߵľ���
		int contain(const Point_3D &point_temp)const;    //�ж��߶ΰ���һ����

	private:
				
		//���ݳ�Ա
        double xm, yn, zl ;       //ֱ�߷��̵ı���ϵ�� (x-x0)/xm=(y-y0)/yn=(z-z0)/zl
};
//---------------------------------------------------------------------------
//������ά����
class Plane_3D
{
	public:

		//���ݳ�Ա
		double coef[4];			//��ʾ�ռ�ƽ�淽�̵��ĸ�ϵ��ax+by+cz+d=0
		bool virtual_plane;		//��ʾƽ���Ƿ��϶�Ϊ��ƽ��(false:��ƽ�棨��������Ϊ(0,0,0)��; true:����ƽ��)

		//���캯��
		Plane_3D(){};
		Plane_3D(double para[]);
		Plane_3D(vector<double> &para);
		Plane_3D(const vector<double> &para);
		
		//��Ա����
        double distance_point_to_plane(const Point_3D *point_temp)const;    //�㵽�ռ�ƽ��ľ���
        double distance_point_to_plane(const Point_3D &point_temp)const;    //�㵽�ռ�ƽ��ľ���
        double distance_point_to_line(const double dx, const double dy, const double dz)const;  //�㵽�ռ�ƽ��ľ���
		int contain(const Point_3D &point_temp)const;    //�жϿռ�ƽ�����һ����
		int contain(const double dx, const double dy, const double dz)const;    //�жϿռ�ƽ�����һ����
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
