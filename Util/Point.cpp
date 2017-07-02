/*
 * Point.cpp
 *
 *  Created on: 09-10-2013
 *      Author: miguel
 */

#include "Point.h"

Point::Point() {
	x = 0;
	y = 0;
}
Point::Point(int x, int y) {
	this->x = x;
	this->y = y;
}
Point::~Point() {

}

long int Point::getX() {
	return x;
}
long int Point::getY() {
	return y;
}
void Point::setX(int x) {
	this->x = x;
}
void Point::setY(int y) {
	this->y = y;
}
void Point::setXY(int x, int y) {
	this->x = x;
	this->y = y;
}
Point * Point::sumar(Point * p) {
	return new Point(this->x + p->x, this->y + p->y);
}
Point * Point::restar(Point * p) {
	return new Point(this->x - p->x, this->y - p->y);
}
