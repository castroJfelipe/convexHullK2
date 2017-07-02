//============================================================================
// Name        : ConvexHullK2Tree.cpp
// Author      : Juan Castro, Gilberto Gutiérrez y Miguel Romero
// Version     :
// Copyright   :
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <stdio.h>
#include "Snapshot/Snapshot.h"
#include <BitSequence.h>
#include <stdlib.h>
#include <queue>
#include <vector>
#include <functional>
#include "Util/TimeMesure.h"
#include "Util/HashMap.h"
#include <unistd.h>
#include <time.h>
#include <algorithm>
#include <string.h>

#define BILLION 1000000000L


using namespace cds_static;
using namespace std;

typedef struct point{
    int x;
    int y;
}POINT;

bool menor(Point a, Point b){
    if (a.getY() != b.getY())
    return a.getY() < b.getY();
    return a.getX() < b.getX();
}

Point pivot;
// returns -1 if a -> b -> c forms a counter-clockwise turn,
// +1 for a clockwise turn, 0 if they are collinear
int ccw(Point a, Point b, Point c) {
	long int area = (b.getX() - a.getX()) * (c.getY() - a.getY()) - (b.getY() - a.getY()) * (c.getX() - a.getX());
	if (area > 0)
	return -1;
	else if (area < 0)
	return 1;
return 0;
}
// returns square of Euclidean distance between two points
long int sqrDist(Point a, Point b) {
long int dx = a.getX() - b.getX();
long int dy = a.getY() - b.getY();
return dx * dx + dy * dy;
}
// used for sorting points according to polar order w.r.t the pivot
bool POLAR_ORDER(Point a, Point b) {
int order = ccw(pivot, a, b);
if (order == 0)
return sqrDist(pivot, a) < sqrDist(pivot, b);
return (order == -1);
}
stack<Point> grahamScan(Point *points, int N) {
	stack<Point> hull;
	if (N < 3)
	return hull;
	// find the point having the least y coordinate (pivot),
	// ties are broken in favor of lower x coordinate
	int leastY = 0;
	for (int i = 1; i < N; i++)
	if (menor(points[i], points[leastY]))
	leastY = i;
	// swap the pivot with the first point
	Point temp = points[0];
	points[0] = points[leastY];
	points[leastY] = temp;
	// sort the remaining point according to polar order about the pivot
	pivot = points[0];
	sort(points + 1, points + N, POLAR_ORDER);
	Point aux;
        aux.setXY(-1,-1);
	hull.push(aux);
	hull.push(points[0]);
	hull.push(points[1]);
	hull.push(points[2]);
    int j = 0, i = 0;
    while(j < (N+3)){
        i = i%N;
        Point top = hull.top();
        hull.pop();
        while (ccw(hull.top(), top, points[i])!= -1) {
         	top = hull.top();
           	hull.pop();
           	if(top.getX() == -1) break;
        }
        hull.push(top);
        hull.push(points[i]);
        j++; i++;
    }
        /*for (int i = 3; i < N; i++) {

	}*/
	return hull;
}


timespec diff(timespec start, timespec end){
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

void printBitString2(BitSequence * a) {
    for (size_t i = 0; i < a->getLength(); i++) {
        fprintf(stdout, "%d", a->access(i));
        if ((i + 1) % 100 == 0) {
            fprintf(stdout, "\n");
        } else if ((i + 1) % 4 == 0) {
            fprintf(stdout, " ");
        }
    }
}

void pruConstruccion(char *pa1, char *pa2, char *pa3, char *pa4) {
//void pruConstruccion() {

    Snapshot * snap;
    //calculo de nivel suponiendo un total de numXYd filas y columnas.

    char a3[200] = {"/home/juan/compacta/workspace/convexHull/convexPuro/src/datos/"};
    strcat(a3,pa2);
    strcat(a3,"/");
    strcat(a3,pa3);
    strcat(a3,"/");
    strcat(a3,pa4);
    strcat(a3,"/");
    strcat(a3,pa1);
    FILE *fp = fopen(a3, "r");

    //FILE *fp = fopen("/home/juan/compacta/workspace/convexHull/convexPuro/src/kCoordenadas.txt", "r");

    int numXYd=16;
    uint edge, objs;
    fscanf(fp,"%d,%d,%d",&numXYd, &edge, &objs);
    cout << "NODOS: " << numXYd << endl;
    uint nodos = numXYd*numXYd;
    uint maxLevel = ceil(log(numXYd) / log(K));
    lkt * tmp = createLKTree(maxLevel);
    uint x,y;
    POINT *pp = new POINT[edge];

    //FILE * f2 = fopen("/home/juan/compacta/workspace/convexHull/convexPuro/registroCNRTamano.txt","a+");
    //fprintf(f2,"Coordenadas: \n");
    for(int j = 0; j < edge; j++){
        fscanf(fp,"%ld,%ld",&y,&x); //Se leen al revés ya que se deben ingresar como si fueran (y,x)
        //fprintf(f2,"%d,%d\n",y,x);
        pp[j].x = x; pp[j].y = y;
        insertNode(tmp,x,y,j);
    }
    uint nx, ny, sx, sy, ex, ey, ox, oy;
    fscanf(fp, "%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld",&nx,&ny,&sx,&sy,&ex,&ey,&ox,&oy);
    //fprintf(f2,"N:(%d,%d)S:(%d,%d)E:(%d,%d)O:(%d,%d)\n",nx,ny,sx,sy,ex,ey,ox,oy);
    snap = createSnapshot(tmp, numXYd, edge, objs);
    fclose(fp);

    char a31[200] = {"/home/juan/compacta/workspace/convexHull/convexPuro/src/HOYTama"};
    strcat(a31,pa2);
    strcat(a31,pa3);
    strcat(a31,pa4);
    strcat(a31,".txt");
    FILE *f2 = fopen(a31,"a+");

    //luego de construir la representación compacta se destruye la temporal!
    destroyLKTree(tmp);
    uint b = 0;
    BitSequence * a = snap->ktree->bt;
    a = snap->ktree->bn;
    a = snap->ktree->bl;
    Point r, aux;
    uint o;
    stack<Point> ch;
    timespec time1, time2;
    Cronometer* cronom = cCronometer();
    //start_clock(cronom);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    //ch = convexHull(snap->ktree, Point(nx,ny), Point(sx,sy), Point(ex,ey), Point(ox,oy), maxLevel);
    int tama = convexHull(snap->ktree, Point(nx,ny), Point(sx,sy), Point(ex,ey), Point(ox,oy), maxLevel);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    double sys_elapsed;
    sys_elapsed = time2.tv_sec + time2.tv_nsec/1E9;
    sys_elapsed -= time1.tv_sec + time1.tv_nsec/1E9;
    //fprintf(stdout,"TIEMPO CHK2: %f\t",sys_elapsed);
    //fprintf(f2,"%f\t",sys_elapsed);

    /*while(!ch.empty()){
    	cout << ch.top().getX() << "\t" << ch.top().getY() << endl;
    	ch.pop();
    }*/

    //fprintf(f2,"%d\t",sizeSnapshot(snap));
    //fprintf(f2,"%d\t",tama);

    fprintf(stdout,"k2tree: %d\t",sizeSnapshot(snap));
    fprintf(stdout,"es. ad: %d\n",tama);

    //fprintf(stdout,"%d\n",tama);

    //double cpu_time = stop_clock(cronom);
    //fprintf(stdout, "TIEMPO DE ALGORITMO CHK2: %f \n",cpu_time);
    //fprintf(f2, "%f \t",cpu_time);

    fclose(f2);

    /*fprintf(stdout, "uint: %d\n",sizeof(uint));
    fprintf(stdout, "int: %d\n",sizeof(int));
    fprintf(stdout, "size_t: %d\n",sizeof(size_t));
    fprintf(stdout, "bitsequence: %d\n",sizeof(BitSequence *));
    fprintf(stdout, "long long int %d\n",sizeof(unsigned long long int *));*/

    //cout << "TAMANO CHK2: " << sizeSnapshot(snap) << endl;
    destroySnapshot(snap);
    delete(pp);
    snap = NULL;
}

void malditoGraham(char *pa1, char *pa2, char *pa3, char *pa4){
	char a2[200] = {"/home/juan/compacta/workspace/convexHull/convexPuro/src/datos/"};
	strcat(a2,pa2);
	strcat(a2,"/");
	strcat(a2,pa3);
	strcat(a2,"/");
	strcat(a2,pa4);
	strcat(a2,"/");
	strcat(a2,pa1);
	FILE *fp = fopen(a2, "r");
	//FILE *fp = fopen("/home/juan/compacta/workspace/convexHull/convexPuro/src/kCoordenadas.txt", "r");
	if(fp == NULL) cout << "Archivo vacio" << endl;
	else cout << "Archivo correcto" << endl;
	int N, auxN;
	fscanf(fp, "%d,%d,%d", &auxN, &N, &auxN);
	Point *points = new Point[N];
        int q1, q2;
        for(int i = 0; i<N; i++){
	  	fscanf(fp,"%d,%d",&q1, &q2);
                points[i].setX(q1);
                points[i].setY(q2);
                if(points[i].getX() == 9 && points[i].getY() == 0) cout << "Indice " << i << endl;
	}
	fclose(fp);
	char a21[200] = {"/home/juan/compacta/workspace/convexHull/convexPuro/src/tiempo"};
	strcat(a21,pa2);
	    strcat(a21,pa3);
	    strcat(a21,pa4);
	    strcat(a21,".txt");
	    FILE *f2 = fopen(a21,"a+");
        /*Cronometer* cronom = cCronometer();
        start_clock(cronom);
        stack<Point> hull = grahamScan(points, N);
        double cpu_time = stop_clock(cronom);
        fprintf(stdout, "TIEMPO DE ALGORITMO MALDITO GRAHAM: %f \n",cpu_time);
        fprintf(f2, "%f \t",cpu_time);*/
	    timespec time1, time2;
	    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
	    stack<Point> hull = grahamScan(points, N);
	    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	    while(!hull.empty()){
	       	cout << hull.top().getX() << "," << hull.top().getY() << endl;
	      	hull.pop();
	    }
	    double sys_elapsed;
	    sys_elapsed = time2.tv_sec + time2.tv_nsec/1E9;
	    sys_elapsed -= time1.tv_sec + time1.tv_nsec/1E9;
	    //fprintf(stdout,"TIEMPO GRAHAM: %f\n",sys_elapsed);
	    fprintf(f2,"%f\t",sys_elapsed);
        //fprintf(f2,"%d\t",sizeof(points)*N);
        //fprintf(f2,"0\n");
        fclose(f2);
        //cout << sizeof(points)*N << endl;
        points = NULL;
        delete points;
}

void sacaPuntos(char *pa1, char *pa2, char *pa3, char *pa4) {
//void sacaPuntos() {

	char a1[200] = {"/home/juan/compacta/workspace/convexHull/convexPuro/src/datos/"};
	strcat(a1,pa2);
	strcat(a1,"/");
	strcat(a1,pa3);
	strcat(a1,"/");
	strcat(a1,pa4);
	strcat(a1,"/");
	strcat(a1,pa1);
	FILE *fp = fopen(a1, "r");

	Snapshot * snap;

    //FILE *fp = fopen("/home/juan/compacta/workspace/convexHull/convexPuro/src/kCoordenadas.txt", "r");

    int numXYd=16;
    uint edge, objs;
    fscanf(fp,"%d,%d,%d",&numXYd, &edge, &objs);
    cout << "NODOS: " << numXYd << endl;
    uint nodos = numXYd*numXYd;
    uint maxLevel = ceil(log(numXYd) / log(K));
    lkt * tmp = createLKTree(maxLevel);
    uint x,y;
    POINT *pp = new POINT[edge];
    for(int j = 0; j < edge; j++){
        fscanf(fp,"%d,%d",&y,&x); //Se leen al revés ya que se deben ingresar como si fueran (y,x)
        //fprintf(f2,"%d,%d\n",y,x);
        pp[j].x = x; pp[j].y = y;
        insertNode(tmp,x,y,j);
    }
    uint nx, ny, sx, sy, ex, ey, ox, oy;
    fscanf(fp, "%d,%d,%d,%d,%d,%d,%d,%d",&nx,&ny,&sx,&sy,&ex,&ey,&ox,&oy);
    //fprintf(f2,"N:(%d,%d)S:(%d,%d)E:(%d,%d)O:(%d,%d)\n",nx,ny,sx,sy,ex,ey,ox,oy);
    snap = createSnapshot(tmp, numXYd, edge, objs);
    fclose(fp);
    //luego de construir la representación compacta se destruye la temporal!
    destroyLKTree(tmp);
    uint b = 0;
    BitSequence * a = snap->ktree->bt;
    a = snap->ktree->bn;
    a = snap->ktree->bl;
    Point r;
    Point * points = new Point[edge];
    uint o;

    stack<Point> ch;
    for (int xa = 0; xa < edge; xa++) {
        points[xa].setXY(-1,-1);
    }

    /*int ni = 0;
    for(ni = 0; ni < edge; ni++){
       	printf("%d,%d\n",points[ni].getX(), points[ni].getY());
    }*/
    points = desintegra(snap->ktree, edge);


    int ka;
    for(ka = 0; ka < edge; ka++){
        if(points[ka].getX() == -1) break;
    }
    /*int na;
    for(na = 0; na < ka; na++){
    	printf("%d,%d\n",points[na].getX(), points[na].getY());
    }*/
    points = NULL;
    delete points;
    Point * po = new Point[ka];

    //FILE *ff = fopen("/home/juan/compacta/workspace/convexHull/convexPuro/src/desintegra.txt","w");


    //FILE * f2 = fopen("/home/juan/compacta/workspace/convexHull/convexPuro/registroCNRTamano.txt","a+");

    char a11[200] = {"/home/juan/compacta/workspace/convexHull/convexPuro/src/HOYTama"};
    strcat(a11,pa2);
    strcat(a11,pa3);
    strcat(a11,pa4);
    strcat(a11,".txt");
	FILE *f2 = fopen(a11,"a+");

	/*Cronometer* cronom = cCronometer();
    start_clock(cronom);
    po = desintegra(snap->ktree,ka);
    stack<Point> hull = grahamScan(po, edge);
    double cpu_time = stop_clock(cronom);

    fprintf(stdout, "TIEMPO DE ALGORITMO SACAPUNTOS: %f \n",cpu_time);
    fprintf(f2, "%f \n",cpu_time);*/

	timespec time1, time2;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
	po = desintegra(snap->ktree,ka+1);
	po[ka].setX(po[ka-1].getX());
	po[ka].setY(po[ka-1].getY());
	stack<Point> hull = grahamScan(po, ka+1);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    double sys_elapsed;
    sys_elapsed = time2.tv_sec + time2.tv_nsec/1E9;
    sys_elapsed -= time1.tv_sec + time1.tv_nsec/1E9;
    //fprintf(stdout,"TIEMPO SACAPUNTOS: %f\n",sys_elapsed);

    //fprintf(f2,"%f\n",sys_elapsed);
    //fprintf(f2,"%d\t",sizeSnapshot(snap));
    //fprintf(stdout,"%d\n",sizeof(po)*ka);

    //fprintf(fp,"%d\n",sizeof(Point)*ka);

    //fprintf(stdout,"%d\n",sizeof(Point));
    //fprintf(stdout,"%d\n",ka);
    fprintf(stdout,"sacapunto: %d\n",sizeof(Point)*ka);
    //fclose(f2);

    //FILE *co = fopen("/home/juan/compacta/workspace/convexHull/convexPuro/src/aa.txt","w");
    /*while(!hull.empty()){
      	//fprintf(co,"%d\t%d\n",hull.top().getX(),hull.top().getY());
    	cout << hull.top().getX() << "\t" << hull.top().getY() << endl;
       	hull.pop();
    }*/
    //fclose(co);

    /*fprintf(f2,"%d\t",sizeSnapshot(snap));
    fprintf(f2,"%d\t",sizeof(po)*ka);
    fclose(f2);*/


    //cout << "TAMANO SACAPUNTOS: " << sizeof(po)*ka << endl;
    po = NULL;
    delete po;
    destroySnapshot(snap);
    delete(pp);
    snap = NULL;
}

void generaAleatorio(){
    srand(time(NULL));
        int coo = 0;
        cout << "Cuantas coordenadas generar?: ";
        cin >> coo;
        int n = 0;
        cout << "Defina el N: ";
        cin >> n;

        POINT p[coo];
        POINT aux;
        POINT norte, sur, este, oeste;
        norte.x = n-1; norte.y = n-1;       sur.x = 0; sur.y = 0;
        este.x = 0; este.y = 0;
        oeste.x = n-1; oeste.y = n-1;
        uint flag = 0;
        FILE *fp = fopen("/home/juan/compacta/workspace/convexHull/convexPuro/src/kCoordenadas.txt", "w");
        fprintf(fp, "%d,%d,%d\n",n,coo,coo+1);
        for(int i = 0; i < coo; i++){
            cout << "Loop " << i << " - Coo: " << coo <<endl;
            aux.x = rand()%(n);
            aux.y = rand()%(n);
            for(int j = 0; j < coo; j++){
                if(aux.x == p[j].x && aux.y == p[j].y){
                    flag = 1;
                    break;
                }
            }
            if(flag == 0){
                p[i] = aux;
                fprintf(fp, "%d,%d\n",aux.x,aux.y);
                if(aux.y >= sur.y) sur = aux;
                if(aux.y <= norte.y) norte = aux;
                if(aux.x >= este.x) este = aux;
                if(aux.x <= oeste.x) oeste = aux;
                cout << "Punto: " << p[i].x << "," << p[i].y << endl;
            } else{
                flag = 0;
                i--;
            }
        }
        fprintf(fp, "%d,%d,%d,%d,%d,%d,%d,%d",norte.x,norte.y,sur.x,sur.y,este.x,este.y,oeste.x,oeste.y);
        fclose(fp);
        printf("Extremos N(%d,%d), S(%d,%d), E(%d,%d), O(%d,%d)",norte.x,norte.y,sur.x,sur.y,este.x,este.y,oeste.x,oeste.y);
}

int main(int argc, char *argv[]) {
	pruConstruccion(argv[4], argv[1], argv[2], argv[3]);
    //malditoGraham(argv[4], argv[1], argv[2], argv[3]);
    //sacaPuntos(argv[4], argv[1], argv[2], argv[3]);
    return 0;
}

/*int main() {
	generaAleatorio();
	pruConstruccion();
    //malditoGraham(argv[4], argv[1], argv[2], argv[3]);
    //sacaPuntos();
	//sacaPuntos();
    return 0;
}*/
