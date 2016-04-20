//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Geometry_2D.cpp
//OBJECTIVE:	The definitions of point, line and shape in 2D
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#include "Geometry_2D.h"

//---------------------------------------------------------------------------
//��ά����
Point_2D::Point_2D( double px, double py )
{
	x = px;
	y = py;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator+( Point_2D &pt )
{
	Point_2D rp( x + pt.x, y + pt.y);
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator+( const Point_2D &pt )
{
	Point_2D rp( x + pt.x, y + pt.y);
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator+( double d )
{
	Point_2D rp( x + d, y + d);
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator-( Point_2D &pt )
{
	Point_2D rp( x - pt.x, y - pt.y );
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator-( const Point_2D &pt )
{
	Point_2D rp( x - pt.x, y - pt.y );
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator-( double d )
{
	Point_2D rp( x - d, y - d );
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator*( double d )
{
	Point_2D rp( x*d, y*d );
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator/( double d )
{
	Point_2D rp( x/d, y/d );
	return rp;
}
//---------------------------------------------------------------------------
bool Point_2D::operator==( Point_2D &pt )
{
	return (x==pt.x&&y==pt.y);
}
//---------------------------------------------------------------------------
bool Point_2D::operator!=( Point_2D &pt )
{
	return (x!=pt.x||y!=pt.y);
}
//---------------------------------------------------------------------------
double Point_2D::distance_to(const Point_2D &pt )const
{
	double rv2 = (x-pt.x)*(x-pt.x)+(y-pt.y)*(y-pt.y);
	return sqrt(rv2);
}
//---------------------------------------------------------------------------
double Point_2D::distance_to(const double &px, const double &py)const
{
	double rv2 = (x-px)*(x-px)+(y-py)*(y-py);
	return sqrt(rv2);
}
//===========================================================================

//��ά����
//---------------------------------------------------------------------------
//���캯��
Line_2D::Line_2D(Point_2D p0, Point_2D p1)
{
	point[0] = p0;
	point[1] = p1;
	A = p1.y - p0.y ;
	B = -(p1.x - p0.x);
	D = -(A*p0.x+B*p0.y);
	len = length();
	if(len==0) virtual_line = false;
	else virtual_line = true;
}
//---------------------------------------------------------------------------   
//�߶εĳ���
double Line_2D::length() 
{
	double dx = point[1].x-point[0].x;
	double dy = point[1].y-point[0].y;
	return sqrt(dx*dx+dy*dy);
}
//---------------------------------------------------------------------------
//�㵽ֱ�ߵľ���
double Line_2D::distance_point_to_line(const Point_2D *point_temp)const
{
	if(A==0&&B==0) hout << "ע�⣡���߶��˻���һ���㣡��㵽ֱ�ߵľ���ʱ����ĸΪ�㣡" <<endl;
	return fabs(A*point_temp->x+B*point_temp->y+D)/sqrt(A*A+B*B);
}
double Line_2D::distance_point_to_line(const Point_2D &point_temp)const
{
	if(A==0&&B==0) hout << "ע�⣡���߶��˻���һ���㣡��㵽ֱ�ߵľ���ʱ����ĸΪ�㣡" <<endl;
	return fabs(A*point_temp.x+B*point_temp.y+D)/sqrt(A*A+B*B);
}
double Line_2D::distance_point_to_line(double dx, double dy)const
{
	if(A==0&&B==0) hout << "ע�⣡���߶��˻���һ���㣡��㵽ֱ�ߵľ���ʱ����ĸΪ�㣡" <<endl;
	return fabs(A*dx+B*dy+D)/sqrt(A*A+B*B);
}
//---------------------------------------------------------------------------
//�㵽ֱ�ߵ�·�����жϵ�����ֱ�ߵ�һ�໹����һ��
double Line_2D::path_point_to_line(const Point_2D &point_temp)const
{
	if(A==0&&B==0) hout << "ע�⣡���߶��˻���һ���㣡��㵽ֱ�ߵľ���ʱ����ĸΪ�㣡" <<endl;
	return A*point_temp.x+B*point_temp.y+D;
}
//---------------------------------------------------------------------------
//�ж��߶ΰ���һ����
int Line_2D::contain(const Point_2D &point_temp)const
{
	//�㵽�߶������˵�ľ�����������˵�֮��ľ���
	if( fabs(point_temp.distance_to(point[0])+point_temp.distance_to(point[1])-point[0].distance_to(point[1]))>Zero ) return 0; 
	return 1;
}
int Line_2D::contain(const double &px, const double &py)const
{
	Point_2D point_temp(px, py);
	//�㵽�߶������˵�ľ�����������˵�֮��ľ���
	if( fabs(point_temp.distance_to(point[0])+point_temp.distance_to(point[1])-point[0].distance_to(point[1]))>Zero ) return 0; 
	return 1;
}
//===========================================================================

//������
//---------------------------------------------------------------------------
//���캯��
Rectangle::Rectangle(Point_2D p0, Point_2D p1, Point_2D p2, Point_2D p3)
{
	point[0] = p0;
	point[1] = p1;
	point[2] = p2;
	point[3] = p3;
	length = p2.x-p0.x;
	width = p2.y-p0.y;
	area = length*width;
	if(length>0&&width>0) virtual_rect = 1;
	else virtual_rect = 0;
}
Rectangle::Rectangle(Point_2D p0, double len, double wid)
{
	point[0].x = p0.x;		point[0].y = p0.y;
	point[1].x = p0.x;		point[1].y = p0.y+wid;
	point[2].x = p0.x+len;		point[2].y = p0.y+wid;
	point[3].x = p0.x+len;		point[3].y = p0.y;
	length = len;
	width = wid;
	area = length*width;
	if(length>0&&width>0) virtual_rect = 1;
	else virtual_rect = 0;
}
Rectangle::Rectangle(Point_2D p0, Point_2D p2)
{
	point[0] = p0;
	point[1].x = p0.x;
	point[1].y = p2.y;
	point[2] = p2;
	point[3].x = p2.x;
	point[3].y = p0.y;
	length = p2.x-p0.x;
	width = p2.y-p0.y;
	area = length*width;
	if(length>0&&width>0) virtual_rect = 1;
	else virtual_rect = 0;
}
//---------------------------------------------------------------------------
//������ε����
double Rectangle::calculate_area()
{
	return length*width;
}
//---------------------------------------------------------------------------
//�жϾ��������Ƿ������һ���㣨����ֵ��0û������1������
int Rectangle::contain(const Point_2D &poi)const
{
	if(  poi.x>=point[0].x && poi.x<=point[2].x &&
		  poi.y>=point[0].y && poi.y<=point[2].y ) return 1;
	else return 0;
}
//---------------------------------------------------------------------------
//�ж�һ�����Ƿ񱻰����ھ��������ڲ�(�߽����)������ֵ��0û������1������
int Rectangle::contain_in(const Point_2D &poi)const
{
	if(  poi.x>point[0].x && poi.x<point[2].x &&
		  poi.y>point[0].y && poi.y<point[2].y ) return 1;
	else return 0;
}
//---------------------------------------------------------------------------
//�жϾ��������Ƿ������һ���߶Σ�����ֵ��0û������1������
int Rectangle::contain(const Line_2D &line)const
{
	for(int i=0; i<2; i++)
	if( contain(line.point[i])==0 ) return 0;
	return 1;
}
//---------------------------------------------------------------------------
//�жϾ��������Ƿ��������һ���Σ�����ֵ��0û������1������
int Rectangle::contain(const Rectangle &rect)const
{
	for(int i=0; i<4; i++)
	if( contain(rect.point[i])==0 ) return 0;
	return 1;
}
//---------------------------------------------------------------------------
//�ж�һ���߶��Ƿ񱻰����ھ����ڲ����ھ��α߽��ϲ��㣨����ֵ��0û������1������ 
int Rectangle::contain_in(const Line_2D &line)const
{
	if( contain(line.point[0])==1&&contain(line.point[1])==1 ) 
	{
		if( line.point[0].x == line.point[1].x&&(line.point[0].x==point[0].x||line.point[0].x==point[2].x) ) return 0;
		if( line.point[0].y == line.point[1].y&&(line.point[0].y==point[0].y||line.point[0].y==point[2].y) ) return 0;
		return 1;
	}
	return 0;
}
//---------------------------------------------------------------------------
//�ж��߶��Ƿ�;����ཻ���ཻ�����߶ε������˵㶼�ھ��εı߽��ϣ�������һ���˵��ھ�����
//���������һ���˵��ھ��α߽��ϣ���һ�˵��ھ����ⲿ�����ཻ)(����ֵ��0���ཻ��1�ཻ)
int Rectangle::overlap(const Line_2D &line)const
{
	if( contain(line.point[0])==1||contain(line.point[1])==1 )
	{
		if( line.point[0].x == line.point[1].x) 
		{
			if( line.point[0].y > line.point[1].y&&(line.point[0].y <= point[0].y||line.point[1].y>=point[2].y) ) return 0;
			if( line.point[0].y < line.point[1].y&&(line.point[1].y <= point[0].y||line.point[0].y>=point[2].y) ) return 0;
		}
		if( line.point[0].y == line.point[1].y) 
		{
			if( line.point[0].x > line.point[1].x&&(line.point[0].x <= point[0].x||line.point[1].x>=point[2].x) ) return 0;
			if( line.point[0].x < line.point[1].x&&(line.point[1].x <= point[0].x||line.point[0].x>=point[2].x) ) return 0;
		}
		return 1;
	}
	return 0;
}
//---------------------------------------------------------------------------
//�������Rect1��Rect2���ص�������
int Rectangle::rectangles_overlapping(const Rectangle &Rect1, const Rectangle &Rect2)
{
	//�������εĹ������֣����������ε���͵�����ߵ�������������ѡ�м������ֵ
	//�������ε���͵���ȡ�������ֵ
	if(Rect1.point[0].x<=Rect2.point[0].x)	point[0].x = Rect2.point[0].x;
	else point[0].x = Rect1.point[0].x;
	if(Rect1.point[0].y<=Rect2.point[0].y)	point[0].y = Rect2.point[0].y;
	else point[0].y = Rect1.point[0].y;
	//�������ε���ߵ���ȡ������Сֵ
	if(Rect1.point[2].x<=Rect2.point[2].x)	point[2].x = Rect1.point[2].x;
	else point[2].x = Rect2.point[2].x;
	if(Rect1.point[2].y<=Rect2.point[2].y)	point[2].y = Rect1.point[2].y;
	else point[2].y = Rect2.point[2].y;
	
	point[1].x = point[0].x;
	point[1].y = point[2].y;
	point[3].x = point[2].x;
	point[3].y = point[0].y;
	length = point[2].x-point[0].x;
	width = point[2].y-point[0].y;
	area = length*width;
	if(length>0&&width>0) virtual_rect = 1;
	else virtual_rect = 0;  //û���ص����Σ����߽����ص�һ���һ���߶�

	return virtual_rect;
}
//---------------------------------------------------------------------------
//���������ڲ�С����Ϊ�������Ź���ָ�
int Rectangle::make_nine_grids_by(const Rectangle rect, vector<Rectangle> &grids)
{
	if( contain(rect)==0 ) 
	{
		hout << "�þ��β��ڽ�Ҫ���ָ�Ĵ���������ڣ� ���飡" << endl;
		return 0;
	}	
	//��¼�����ķָ��
	Point_2D point_base[4];
	point_base[0] = point[0];
	point_base[1] = rect.point[0];
	point_base[2] = rect.point[2];
	point_base[3] = point[2];
	//����Ź��������еĽڵ㣨�������ң��������ϣ�
	Point_2D point_ex[16];
	for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
		{
			point_ex[i*4+j].x = point_base[j].x;
			point_ex[i*4+j].y = point_base[i].y;
		}
	//��¼����ĸ����Ӿ���
	grids.clear();  //�������
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
		{
			Rectangle rect_temp( point_ex[i*4+j], point_ex[(i+1)*4+j], point_ex[(i+1)*4+j+1], point_ex[i*4+j+1] );
			grids.push_back(rect_temp);
		}

	return 1;
}

//---------------------------------------------------------------------------
//������ε�ȫ����Ϣ�������ĸ��������꣬���������������α�ǣ�
void Rectangle::output_parameters()
{
	hout << "The coordinates of four corner points (starting from lower-left corner point in clockwise direction):" << endl;
	for(int i=0; i<4; i++) hout << "(" << point[i].x << ", " << point[i].y << ")" << "   ";
	hout << endl;
	hout << "The length, width, area and virtual_rect tag:" << endl;
	hout << "Length " << length << ", Width " << width << ", area " << area << ", virtual_rect " << virtual_rect << endl;
}
//===========================================================================
