/*
 * Point.h
 *
 *  Created on: 09-10-2013
 *      Author: miguel
 */

#ifndef POINT_H_
#define POINT_H_

class Point {
public:
	Point();
	Point(int x, int y);
	long int getX();
	long int getY();
	void setX(int x);
	void setY(int y);
	void setXY(int x, int y);
	Point * sumar(Point * p);
	Point * restar(Point * p);

	virtual ~Point();
private:
	long int x;
	long int y;
};

#endif /* POINT_H_ */
