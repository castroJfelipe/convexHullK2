/*
 * LKTree.h
 *  Es una extensión de los ktree que permite
 *  el etiquetado de aristas.
 *  la idea es poder etiquetar con una o más etiquetas una arista.
 *  cada valor de 32 bits.
 *  Created on: 15-11-2012
 *      Author: Miguel Romero Vásquez
 */

#ifndef LKTREE_H_
#define LKTREE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stack>

#include <BitSequence.h>
using namespace cds_static;

#include "../Util/LinkedList.h"
#include "../Util/ListIterator.h"
#include "../Util/LinkedList.h"
#include "../Util/HashMap.h"
#include "../Util/Queue.h"
#include "../Util/ObjectTimePosition.h"
#include "../Util/Point.h"
#include "LabelsMap.h"

#ifndef uint
#define uint unsigned int
#endif

#ifndef K
#define K 2
#endif

void detectMemoryCorrupt();

//Representación temporal usada en la construcción del arbol
//pero no es compacta.

typedef struct lnode {
    uint data;
    //lista de etiquetas
    LinkedList ** labels;
    struct lnode** child;
} L_NODE;

typedef struct slkt {
    L_NODE * root;
    int max_Level;
    uint numberNodes;
    uint numberLeaves;
    uint numberTotalLeaves;
    uint totalLabels;
} lkt;
L_NODE * createL_Node();
lkt * createLKTree(uint maxlevels);
void insertNode(lkt * tree, int x, int y, uint label);
void _insertNode(lkt * tree, int x, int y, uint * labelArray);
int destroyLKTree(lkt * tree);


typedef struct matrixRep2
{
    BitSequence * bt;       //Bitmap representando el arbol
    BitSequence * bn;       //Bitmap representando el padre de las hojas
    BitSequence * bl;       //Bitmap representando las hojas
    int maxLevel;           //Nivel maximo del arbol
    size_t numberOfNodes;
    size_t numberOfEdges;
    unsigned long long int  * div_level_table;
}MREP2;

//reprecentación compacta de un lkt.
typedef struct sSnapshot {
    uint totObj;
    MREP2 * ktree;
    LabelsMap * labels;
} Snapshot;

typedef struct puntosD{
	double x;
	double y;
}PUNTOSD;

Snapshot * createSnapshot(lkt * tree, uint numberOfNodes,
                                uint numberOfEdges, uint totObject);



LinkedList *  rangeQuery(Snapshot * rep, uint p1, uint p2, uint q1, uint q2);
void rangeQuery(Snapshot * snap, uint p1, uint p2, uint q1, uint q2, uint * &Oid, uint * &X, uint * &Y, uint &n);
Point * getObjectPos(Snapshot * rep, uint oid);
void NN(Snapshot * rep, Point  q, uint &oid, Point & p);
void kNN(Snapshot * rep, Point q, uint & kn, uint * oid, Point * p);
void kNNFilter(Snapshot * rep, Point q, uint maxDesp ,uint & kn, uint * oid, Point * p);
void destroySnapshot(Snapshot * rep);

/**
 * convex hull
 *      Encuentra la cerradura convexa de los puntos que existen en el k2tree.
 * Parámetros:
 * Entrada:
 *       k2tree: representación compacta del espacio con un k^2-Tree
 *       n: punto más al norte
 *       s: punto más al sur
 *       e: punto más al este
 *       o: punto más al oeste
 * Salida:
 *      ch: arreglo con todos los puntos que forma el convex hull.
 */


//stack<Point> convexHull(MREP2 * k2tree, Point n, Point s, Point e, Point o, uint maxlevel);
int convexHull(MREP2 * k2tree, Point n, Point s, Point e, Point o, uint maxlevel);
Point * desintegra(MREP2 * k2tree, int n);


size_t sizeSnapshot(Snapshot * lktrep);
size_t sizeMREP(MREP2 * rep);
#endif /* LKTREE_H_ */
