/*
 * LKTree.cpp
 *
 *  Created on: 15-11-2012
 *      Author: miguel
 */

#include "Snapshot.h"
#include "../Util/HashIter.h"
#include "../Util/Factory.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <stack>
#include <queue>
#include <algorithm>
#include <functional>
#include <time.h>

using namespace std;

#define MAX_RESULT_LIMIT 200000
//void detectMemoryCorrupt() {
//  int * i = (int*) malloc(sizeof(int));
//  free(i);
//  int * j = (int*) malloc(sizeof(int));
//  free(j);
//  int * k = (int*) malloc(sizeof(int));
//  free(k);
//  int * l = (int*) malloc(sizeof(int));
//  free(l);
//}

//Operacioens con representación de árbol
//creación de la representación
L_NODE * createL_Node() {
    L_NODE * resp = (L_NODE *) malloc(sizeof(L_NODE));
    resp->data = 0;
    resp->child = NULL;
    resp->labels = NULL;
    return resp;
}

lkt * createLKTree(uint maxlevels) {
    lkt * tree = (lkt *) malloc(sizeof(lkt));
    tree->root = createL_Node();
    tree->max_Level = maxlevels;
    tree->numberNodes = 0;
    tree->numberLeaves = 0;
    tree->numberTotalLeaves = 0;
    tree->totalLabels = 0;
    return tree;
}

//delete a sub-tree recursively
int destroyLNode(L_NODE * node) {
    int resp = 0;
    if (node == NULL) {
        return resp;
    } else {
        if (node->child != NULL) {
            //delete child recursively
            for (int i = 0; i < K * K; i++) {
                resp += destroyLNode(node->child[i]);
                node->child[i] = NULL;
            }
            free(node->child);
            node->child = NULL;
            resp += sizeof(L_NODE*) * K * K + sizeof(L_NODE*) * K * K;
        }
        if (node->labels != NULL) {
            for (int i = 0; i < K * K; i++) {
                if (node->labels[i] != NULL) {
                    resp += sizeof(LinkedList*)
                            + node->labels[i]->size
                                    * (sizeof(Node) + sizeof(Element)
                                            + sizeof(uint));
                    destroyLinkedList(node->labels[i]);
                    free(node->labels[i]);
                    node->labels[i] = NULL;

                }
            }
            free(node->labels);
            node->labels = NULL;
        }
        free(node);
        node = NULL;
        resp += sizeof(L_NODE) + sizeof(L_NODE*);
    }
    return resp;
}

//elimina la reprecentación temporal de un lktree.
//retorna el total de bytes liberados
int destroyLKTree(lkt * tree) {
    int resp = 0;
    if (tree == NULL)
        return resp;
    //Elimino la raiz
    resp = destroyLNode(tree->root);
    tree->root = NULL;
    //elimino el nodo cabecera.
    free(tree);
    tree = NULL;
    return resp;
}

void insertNode(lkt * tree, int x, int y, uint label) {
    uint * pt = (uint*) malloc(sizeof(uint));
    *pt = label;
    _insertNode(tree, x, y, pt);

}

// guardar aquí información como rangoX[max..min], rangoY[max..min],
//rangoLabel[max..min] y usar esta información para ajustar los limites
//del k2tree, en vez de usar parámetros. asi es más simple la estructura
//con menos parametros.
//insertaNodo con la etiqueta respectiva

void _insertNode(lkt * tree, int x, int y, uint * labelArray) {
    uint i, node = 0;
    unsigned long long int div_level;
    int l = 0;
    L_NODE * n = tree->root;
    while (l <= tree->max_Level) {

        div_level = pow(K, tree->max_Level - l);
        node = (x / div_level) * K + y / div_level;
        if (l == tree->max_Level) {
            if (n->data == 0) {
                tree->numberLeaves++;
            }
            n->data = n->data | (0x1 << node);   //1;
        } else {
            if (n->child == NULL) {
                if (l < tree->max_Level - 1) {
                    tree->numberNodes += K * K;
                } else {
                    tree->numberTotalLeaves += K * K;
                }
                n->child = (L_NODE **) calloc(K * K, sizeof(L_NODE *));
                for (i = 0; i < K * K; i++) {
                    n->child[i] = NULL;
                }
            }
            if (n->child[node] == NULL) {
                n->child[node] = createL_Node();
            }

            //n apunta al nodo insertado (intermedio)
            n = n->child[node];
        }
        x = x % div_level;
        y = y % div_level;
        l++;
    }
    //aquí n es el último nodo insertado, por lo tanto la hoja.
    //hay que insertar la etiqueta en la lista asosiada al nodo
    //
    tree->totalLabels++;

    if (n->labels == NULL) {
        n->labels = (LinkedList**) calloc(K * K, sizeof(LinkedList*));
    }
    if (n->labels[node] == NULL) {
        //se crea la lista de etiquetas si esta no existe
        n->labels[node] = cLinkedList();
    }
    //aquí es importante que el key de element sea node, este
    //valor será usado en la construcción de la representación
    //compacta.
    add(n->labels[node], cElementKey(node, labelArray));

}

//Operaciones con representación compacta

Snapshot * createSnapshot(lkt * tree, uint numberOfNodes, uint numberOfEdges,
uint totObject) {
    //NO destruye lkt, hay que hacerlo desde afuera.
    //Estrategia general:
    //usando el arbol se construye el ktree de los puntos, y se construye un
    //arbol para las etiquetas.  Luego se construye el ktree para las estiquetas
    //con dicho arbol.

    //en este caso los valores de Numero de objetos y numero de arcos no se
    //si son necesarios, hay que revisar.

    Snapshot * Lrep;
    Lrep = (Snapshot *) malloc(sizeof(Snapshot));
    Lrep->totObj = totObject;
    Lrep->ktree = (MREP2 *) malloc(sizeof(MREP2));
    Lrep->ktree->maxLevel = tree->max_Level;
    Lrep->ktree->numberOfNodes = numberOfNodes;
    Lrep->ktree->numberOfEdges = numberOfEdges;
    size_t bits_BT_len = tree->numberNodes;
    size_t bits_BN_len = tree->numberTotalLeaves;
    size_t bits_LI_len = tree->numberLeaves * K * K;

    LinkedList *labelLists = cLinkedList();

    uint * bits_BT = new uint[uint_len(bits_BT_len, 1)]();
    uint * bits_BN = new uint[uint_len(bits_BN_len, 1)]();
    uint * bits_LI = new uint[uint_len(bits_LI_len, 1)]();

    uint k, j, queuecont, conttmp, node, pos = 0;

    unsigned long long div_level;
    int i;
    char isroot = 1;
    Queue * q = createEmptyQueue();
    L_NODE * subTree;
    enqueue(q, tree->root);
    queuecont = 1;

    for (i = 0; i < tree->max_Level; i++) {
        //i:iteración por nivel hasta el penúltimo
        conttmp = 0;
        div_level = pow(K, tree->max_Level - i);
        for (k = 0; k < queuecont; k++) {
            //por elemento de la cola de un mismo nivel
            subTree = (L_NODE *) dequeue(q);
            if ((subTree != NULL) && (subTree->child != NULL)) {
                for (j = 0; j < K * K; j++) {
                    node = j;
                    conttmp++;
                    enqueue(q, subTree->child[node]);
                }
                if (!isroot)
                    bitset(bits_BT, pos);
            }
            if (!isroot)
                pos++;
            isroot = 0;
        }
        queuecont = conttmp;
    }

    Lrep->ktree->bt = PlainBitSequenceFactory(bits_BT, bits_BT_len);
    pos = 0;
    size_t pos_inf = 0;

    //  fprintf(stderr,"Empezando bitmap de hojas utiles\n");
    while (!queueIsEmpty(q)) {
        subTree = (L_NODE *) dequeue(q);
        if ((subTree != NULL) && (subTree->data)) {
            //si no es 0
            bitset(bits_BN, pos);
            //hay que separar las listas

            for (i = 0; i < K * K; i++) {
                if ((subTree->data) & (0x1 << i)) {
                    bitset(bits_LI, pos_inf);
                    add(labelLists, cElement(subTree->labels[i]));
                    labelLists->count += subTree->labels[i]->size;
                }
                pos_inf++;
            }
        }
        pos++;
    }
    destroyQueue(q);
    Lrep->ktree->bn = PlainBitSequenceFactory(bits_BN, bits_BN_len);
    Lrep->ktree->bl = PlainBitSequenceFactory(bits_LI, bits_LI_len);

    Lrep->ktree->div_level_table = (long long uint *) malloc(
            sizeof(long long uint) * (Lrep->ktree->maxLevel + 1));
    for (i = 0; i <= Lrep->ktree->maxLevel; i++)
        Lrep->ktree->div_level_table[i] = pow(K, Lrep->ktree->maxLevel - i);

    //compactando las listas de id.
    Lrep->labels = new LabelsMap(labelLists, totObject);
    //Eliminando los nodos de las listas, pero sin eliminar las listas de ID,
    //para evitar un doble free al eliminar el arbol temporal que también
    //apunta a las mismas listas de objetos.
    destroyLinkedListWhithoutElementPointer(labelLists);
    free(labelLists);
    labelLists = NULL;

    delete[] bits_BT;
    bits_BT = NULL;
    delete[] bits_BN;
    bits_BN = NULL;
    delete[] bits_LI;
    bits_LI = NULL;
    return Lrep;

}

void destroySnapshot(Snapshot * rep) {
    //destroyRepresentation(rep->ktree);
    delete rep->ktree->bl;
    delete rep->ktree->bn;
    delete rep->ktree->bt;
    free(rep->ktree->div_level_table);
    rep->ktree->bl = NULL;
    rep->ktree->bn = NULL;
    rep->ktree->bt = NULL;
    rep->ktree->div_level_table = NULL;
    free(rep->ktree);
    rep->ktree = NULL;

    delete rep->labels;
    rep->labels = NULL;
    free(rep);
    rep = NULL;
}

void recursiveGetRank(Snapshot * snap, uint p1, uint p2, uint q1, uint q2,
uint dp,
                        uint dq, int x, int l, LinkedList * rs) {
    MREP2 * rep = snap->ktree;
    uint i = 0, j, leaf;
    uint y, p1new, p2new, q1new, q2new;
    unsigned long int divlevel;
    if (l == rep->maxLevel) {
        //recorrido por el bitarray leavesInf
        leaf = x + i * p1;
        for (i = p1; i <= p2; i++) {
            for (j = q1; j <= q2; j++) {
                leaf = x + j;
                if (rep->bl->access(leaf)) {
                    //dp + i => posición X del objeto econtrado
                    //dq + j => posición Y del objeto encontrado.
                    snap->labels->getObjects(rep->bl->rank1(leaf), dp + i,
                                                dq + j, rs);
                }
            }
            leaf += K;
        }

    }

    if ((l == rep->maxLevel - 1)
            && (rep->bn->access(x - rep->bt->getLength()))) {
        //recorrido por el bitarray bn
        y = (rep->bn->rank1(x - rep->bt->getLength()) - 1) * K * K;
        for (i = p1; i <= p2; i++) {
            for (j = q1; j <= q2; j++) {
                recursiveGetRank(snap, 0, 0, 0, 0, dp + i, dq + j,
                                    y + K * i + j, l + 1, rs);
            }
        }

    }
    if ((x == -1) || ((l < rep->maxLevel - 1) && (rep->bt->access(x)))) {
        //recorrido por el bitarray bt
        y = (x == -1) ? 0 : rep->bt->rank1(x) * K * K;
        divlevel = rep->div_level_table[l + 1];
        for (i = p1 / divlevel; i <= p2 / divlevel; i++) {
            p1new = 0;
            if (i == p1 / divlevel)
                p1new = p1 % divlevel;
            p2new = divlevel - 1;
            if (i == p2 / divlevel)
                p2new = p2 % divlevel;
            for (j = q1 / divlevel; j <= q2 / divlevel; j++) {
                q1new = 0;
                if (j == q1 / divlevel)
                    q1new = q1 % divlevel;
                q2new = divlevel - 1;
                if (j == q2 / divlevel)
                    q2new = q2 % divlevel;
                recursiveGetRank(snap, p1new, p2new, q1new, q2new,
                                    dp + divlevel * i, dq + divlevel * j,
                                    y + K * i + j, l + 1, rs);
            }
        }

    }
}

LinkedList * rangeQuery(Snapshot * snap, uint p1, uint p2, uint q1, uint q2) {
    LinkedList * rs = cLinkedList();
    recursiveGetRank(snap, p1, p2, q1, q2, 0, 0, -1, -1, rs);
    return rs;
}

void recursiveGetRank(Snapshot * snap, uint p1, uint p2, uint q1, uint q2,
uint dp,
                        uint dq, int x, int l, uint * &Oid, uint * &X,
                        uint * &Y,
                        uint &n) {
    MREP2 * rep = snap->ktree;
    uint i = 0, j, leaf;
    uint y, p1new, p2new, q1new, q2new;
    unsigned long int divlevel;
    if (l == rep->maxLevel) {
        //recorrido por el bitarray leavesInf
        leaf = x + i * p1;
        for (i = p1; i <= p2; i++) {
            for (j = q1; j <= q2; j++) {
                leaf = x + j;
                if (rep->bl->access(leaf)) {
                    //dp + i => posición X del objeto econtrado
                    //dq + j => posición Y del objeto encontrado.
                    snap->labels->getObjects(rep->bl->rank1(leaf), dp + i,
                                                dq + j, Oid, X, Y, n);
                }
            }
            leaf += K;
        }

    }

    if ((l == rep->maxLevel - 1)
            && (rep->bn->access(x - rep->bt->getLength()))) {
        //recorrido por el bitarray bn
        y = (rep->bn->rank1(x - rep->bt->getLength()) - 1) * K * K;
        for (i = p1; i <= p2; i++) {
            for (j = q1; j <= q2; j++) {
                recursiveGetRank(snap, 0, 0, 0, 0, dp + i, dq + j,
                                    y + K * i + j, l + 1, Oid, X, Y, n);
            }
        }

    }
    if ((x == -1) || ((l < rep->maxLevel - 1) && (rep->bt->access(x)))) {
        //recorrido por el bitarray bt
        y = (x == -1) ? 0 : rep->bt->rank1(x) * K * K;
        divlevel = rep->div_level_table[l + 1];
        for (i = p1 / divlevel; i <= p2 / divlevel; i++) {
            p1new = 0;
            if (i == p1 / divlevel)
                p1new = p1 % divlevel;
            p2new = divlevel - 1;
            if (i == p2 / divlevel)
                p2new = p2 % divlevel;
            for (j = q1 / divlevel; j <= q2 / divlevel; j++) {
                q1new = 0;
                if (j == q1 / divlevel)
                    q1new = q1 % divlevel;
                q2new = divlevel - 1;
                if (j == q2 / divlevel)
                    q2new = q2 % divlevel;
                recursiveGetRank(snap, p1new, p2new, q1new, q2new,
                                    dp + divlevel * i, dq + divlevel * j,
                                    y + K * i + j, l + 1, Oid, X, Y, n);
            }
        }

    }
}
//arreglo donde la posición 0 tiene el total de elementos (n) y comienza en 1 y
//termina en n. n es un parámetro de salida.
void rangeQuery(Snapshot * snap, uint p1, uint p2, uint q1, uint q2,
uint * &Oid,
                uint * &X, uint * &Y, uint &n) {
    Oid = new uint[MAX_RESULT_LIMIT]();
    X = new uint[MAX_RESULT_LIMIT]();
    Y = new uint[MAX_RESULT_LIMIT]();
    n = 0;
    recursiveGetRank(snap, p1, p2, q1, q2, 0, 0, -1, -1, Oid, X, Y, n);
}

std::stack<uint> * findPath(Snapshot * snap, uint oid) {
    stack<uint> * path = new stack<uint>();
    size_t p = snap->labels->getListID(oid);
    if (p > snap->ktree->bl->countOnes()) {
        //significa que el objeto no está en el snapshot.
        return path;
    }
    size_t x = snap->ktree->bl->select1(p);
    uint i = x % (K * K);                   //posición como hijo en bl
    path->push(i);
    //posicion del padre de x. en bn
    x = snap->ktree->bn->select1((x / (K * K)) + 1);
    x = x + snap->ktree->bt->getLength();
    while (x > 0) {
        i = x % (K * K);                    //posición como hijo
        path->push(i);
        //posicion del padre de x.
        x = (x < (K * K) ? 0 : snap->ktree->bt->select1((x / (K * K))));
    }
    return path;
}

Point* getObjectPos(Snapshot * snap, uint oid) {
    if (oid > snap->totObj) {
        return NULL;
    }
    uint x1 = 0, y1 = 0, x2 = (snap->ktree->numberOfNodes - 1), y2 = x2;
    stack<uint> * path = findPath(snap, oid);
    uint i, colum, row, piv_x, piv_y;
    while (!path->empty()) {
        i = path->top();
        colum = i % K;
        row = i / K;
        piv_x = ((x2 - x1) + 1) / K;
        piv_y = ((y2 - y1) + 1) / K;
        x2 = piv_x * (row + 1) + x1 - 1;
        x1 = x2 - piv_x + 1;
        y2 = piv_y * (colum + 1) + y1 - 1;
        y1 = y2 - piv_y + 1;
        path->pop();
    }
    delete path;
    path = NULL;
    return new Point(x1, y1);
}

typedef struct sArea {
    unsigned long int x1, x2, y1, y2;
} SpatialArea;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * element Queue
 */
typedef struct sElementQueue {
    long int sTree;
    SpatialArea area;
    double priority;
} ElementQueue;

typedef struct chElementQeue{
    long int sTree;
    uint x;
    uint y;
    short int h;
    double distancia;
} CHElementQueue;

typedef struct sCandidateQueue {
    long int oid;
    long int x;
    long int y;
    double priority;
} CandidateElement;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Max Heap comparator
 */
struct maxHeapComparator {
    bool operator()(CandidateElement i, CandidateElement j) {
        return i.priority < j.priority;
    }
};
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Max Heap comparator
 */
struct minHeapComparator {
    bool operator()(ElementQueue i, ElementQueue j) {
        return i.priority > j.priority;
    }
};

double maxDist(Point p, SpatialArea a) {
    long int dx1 = (p.getX() - a.x1);
    dx1 = dx1 >= 0 ? dx1 : -dx1;

    long int dx2 = (p.getX() - a.x2);
    dx2 = dx2 >= 0 ? dx2 : -dx2;

    long int dy1 = p.getY() - a.y1;
    dy1 = dy1 >= 0 ? dy1 : -dy1;

    long int dy2 = p.getY() - a.y2;
    dy2 = dy2 >= 0 ? dy2 : -dy2;

    long int dx = max(dx1, dx2);
    long int dy = max(dy1, dy2);
    return sqrt((dx * dx) + (dy * dy));
}

double minDist(Point p, SpatialArea a) {
    long int x = p.getX();
    long int y = p.getY();
    if (x < a.x1) {
        //hay 3 casos o bien está dentro del rango y1,y2 o fuera.
        if (y < a.y1) {
            long int dx = a.x1 - x;
            long int dy = a.y1 - y;
            return sqrt((dx * dx) + (dy * dy));
        } else if (y <= a.y2) {
            return double(a.x1 - x);
        } else {
            long int dx = a.x1 - x;
            long int dy = y - a.y2;
            return sqrt((dx * dx) + (dy * dy));
        }
    } else if (x <= a.x2) {
        if (y < a.y1) {
            return double(a.y1 - y);
        } else if (y <= a.y2) {
            return min(min((x - a.x1), (a.x2 - x)), min((y - a.y1), (a.y2 - y)));
        } else {
            return double(y - a.y2);
        }
    } else {
        if (y < a.y1) {
            long int dx = x - a.x2;
            long int dy = a.y1 - y;
            return sqrt((dx * dx) + (dy * dy));
        } else if (y <= a.y2) {
            return double(x - a.x2);
        } else {
            long int dx = x - a.x2;
            long int dy = y - a.y2;
            return sqrt((dx * dx) + (dy * dy));
        }
    }
    return 0.0;
}
inline int access(MREP2 * mrep, long int x) {
    long int btLen = mrep->bt->getLength();
    long int bnLen = mrep->bn->getLength();
    if (x < btLen) {
        return mrep->bt->access(x);
    } else if (x < (btLen + bnLen)) {
        return mrep->bn->access(x - btLen);
    } else {
        //es una hoja de último nivel.
        return mrep->bl->access(x - btLen - bnLen);
    }
}
inline int isLeaf(MREP2 * mrep, long int index) {
    //index está mapeado como si la implementación fuera un único bitmap
    //con la finalidad de encapsular las diferentes implementaciones existentes
    //del k^2-tree
    //verdadero si es una hoja de último nivel sino falso.
    return ((mrep->bt->getLength() + mrep->bn->getLength()) <= index);
}
inline int posInLeaf(MREP2 * mrep, long int index) {
    return index - (mrep->bt->getLength() + mrep->bn->getLength());
}
inline long int firstChild(MREP2 * mrep, long int x) {
    long int btLen = mrep->bt->getLength();
    long int bnLen = mrep->bn->getLength();
    if (x == -1) {
        return 0;
    } else if (x < btLen) {
        return (mrep->bt->access(x) == 0) ? -1 : mrep->bt->rank1(x) * K * K;
    } else if (x < (btLen + bnLen)) {
        return (mrep->bn->access(x - btLen) == 0) ?
                -1 : ((mrep->bt->countOnes() + (mrep->bn->rank1(x - btLen))) * K * K);
    } else {
        //es una hoja de último nivel.
        return -1;
    }
}
//Esto está implementada para un K=2, pero hay que hacerlo para cualquier
//valor de K.
inline SpatialArea getSubArea(SpatialArea a, int i) {
    SpatialArea resp;
    switch (i) {
    case 0:
        resp.x1 = a.x1;
        resp.x2 = (a.x1 + a.x2) / 2;
        resp.y1 = a.y1;
        resp.y2 = (a.y1 + a.y2) / 2;
        break;
    case 1:
        resp.x1 = a.x1;
        resp.x2 = (a.x1 + a.x2) / 2;
        resp.y1 = ((a.y1 + a.y2) / 2) + 1;
        resp.y2 = a.y2;

        break;
    case 2:
        resp.x1 = ((a.x1 + a.x2) / 2) + 1;
        resp.x2 = a.x2;
        resp.y1 = a.y1;
        resp.y2 = (a.y1 + a.y2) / 2;

        break;
    case 3:
        resp.x1 = ((a.x1 + a.x2) / 2) + 1;
        resp.x2 = a.x2;
        resp.y1 = ((a.y1 + a.y2) / 2) + 1;
        resp.y2 = a.y2;
        break;
    }
    return resp;

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*
 * oid y p son arreglos no nulos de tamaño kn
 */
void kNN(Snapshot * rep, Point q, uint & kn, uint * oid, Point * p) {
    priority_queue<CandidateElement, std::vector<CandidateElement>,
            maxHeapComparator> candidates;
    priority_queue<ElementQueue, std::vector<ElementQueue>, minHeapComparator> pQueue;
    ElementQueue e;
    CandidateElement ce;
    e.sTree = 0;
    e.area.x1 = 0;
    e.area.x2 = pow(2, rep->ktree->maxLevel) - 1;
    e.area.y1 = e.area.x1;
    e.area.y2 = e.area.x2;
    e.priority = maxDist(q, e.area);
    long int subTreeIndex;
    double minDis2Child;
    long uint lenCandidates = 0;
    uint numObj = 0;
    uint toAdd = 0;
    uint len = 0;
    uint * resp;
    int rankLeaf;
    pQueue.push(e);
    while (!pQueue.empty()) {
        e = pQueue.top();
        pQueue.pop();
        if (isLeaf(rep->ktree, e.sTree)) {
            //es una hoja de 'ultimo nivel
            lenCandidates = candidates.size();
            if (lenCandidates < kn) {
                numObj = kn - lenCandidates;
                rankLeaf = rep->ktree->bl->rank1(
                        posInLeaf(rep->ktree, e.sTree));
                //trata de recuperar el máximo de objetos para llenar la cola
                //de candidatos.  Puede que no exista más de un objeto en la
                //celda, o bien el número sea menor a los espacios disponibles,
                //en este caso numObj despúes de la llamada tendrá
                //el número real de objetos recuperados
                resp = rep->labels->getUpToKObjects(rankLeaf, numObj);
                for (uint i = 0; i < numObj; i++) {
                    ce.oid = resp[i];
                    //si hubiera información adicional de donde está el objeto
                    //dentro de la celda habría que ajustar la distancia a q
                    //y tambien su x,y.
                    //abría que ajustar la medida.
                    ce.priority = e.priority;
                    //al llegar a la hoja x1==x2 e y1==y2.
                    ce.x = e.area.x1;
                    ce.y = e.area.y1;
                    candidates.push(ce);
                }
                //elimino resp
                delete[] resp;
                resp = NULL;
            } else if (e.priority < (candidates.top().priority)) {
                int i = 0;
                numObj = kn;
                rankLeaf = rep->ktree->bl->rank1(
                        posInLeaf(rep->ktree, e.sTree));
                //aquí no puedo saber cuantos objetos a priory son
                //menores que los actuales candidatos por lo tanto
                //tengo que recuperar Kn objetos, porque en el peor
                //los objetos de la celda reemplazarán a todos los
                //candidatos.  No es necesario recuparar más de Kn
                resp = rep->labels->getUpToKObjects(rankLeaf, numObj);
                do {
                    candidates.pop();
                    ce.oid = resp[i++];
                    ce.priority = e.priority;
                    ce.x = e.area.x1;
                    ce.y = e.area.y1;
                    candidates.push(ce);
                } while (i < numObj && e.priority < (candidates.top().priority));
                delete[] resp;
                resp = NULL;
            }
        } else {
            subTreeIndex = firstChild(rep->ktree, e.sTree);
            if (subTreeIndex != -1) {
                //no es una hoja intermedia, es decir una región sin objetos
                for (int i = 0; i < (K * K); i++) {
                    ElementQueue child;
                    child.sTree = subTreeIndex + i;
                    child.area = getSubArea(e.area, i);
                    child.priority = maxDist(q, child.area);
                    minDis2Child = minDist(q, child.area);

                    if (!((candidates.size() >= kn)
                            && (minDis2Child >= candidates.top().priority))) {
                        pQueue.push(child);
                    }
                }
            }
        }
    }
    kn = candidates.size();
    //puede ser que existan menos kn que los originalmente pedidos
    //entonces en la respuesta kn serán los realmente obtenidos.
    //lo siguiente es un overhead que se podría evitar si se retorna la cola
    //directamente.
    uint i = kn - 1;
    while (!candidates.empty()) {
        ce = candidates.top();
        oid[i] = ce.oid;
        p[i] = Point(ce.x, ce.y);
        i--;
        candidates.pop();
    }
}
//-----------------------------------------------
//para depuración:
void printPush(ElementQueue e) {
    cout << "Push: {" << e.sTree << "," << e.priority << ", (" << e.area.x1
            << ", " << e.area.x2 << ", " << e.area.y1 << ", " << e.area.y2
            << ")}" << endl;
}
void printPop(ElementQueue e) {
    cout << "Pop: {" << e.sTree << "," << e.priority << ", (" << e.area.x1
            << ", " << e.area.x2 << ", " << e.area.y1 << ", " << e.area.y2
            << ")}" << endl;
}
//-------------------------------------------------------

//entrada.
// rep: un snapshot
// q: el punto de referencia
//salida:
//  oid: el identificador del punto más cercano
//  p: la posición del punto más cercanco.

void NN(Snapshot * rep, Point q, uint & oid, Point & p) {
    priority_queue<ElementQueue, std::vector<ElementQueue>, minHeapComparator> pQueue;
    ElementQueue e, c;
    //hay que revisar si al construir el k2tree la posición 0 representa
    //toda el área inicial o bien es el primer hijo de la raiz.
    //por implementación el k2tree construye 1000 como una raiz padre del árbol
    e.sTree = 0;
    e.area.x1 = 0;
    //el 2 es porque k = 2, reemplazar por K.
    e.area.x2 = pow(2, rep->ktree->maxLevel) - 1;
    e.area.y1 = e.area.x1;
    e.area.y2 = e.area.x2;
    e.priority = maxDist(q, e.area);
    c = e;
    c.priority = e.priority + 1;   //explicar el + 1 con un comentario!!!!
    int foundNN = 0;
    long int subTreeIndex;
    double minDis2Child;
    pQueue.push(e);
    printPush(e);
    while (!pQueue.empty()) {
        e = pQueue.top();
        pQueue.pop();
        printPop(e);
        if (isLeaf(rep->ktree, e.sTree)) {
            //si es hoja de último nivel
            if ((e.priority < c.priority)) {
                c = e;
                foundNN = 1;
                //Al encontrar un candidato, puede haber otro mejor
                //en otra área, hay que seguir buscando!
                //aquí se podría hacer una poda del heap tengo que revisar
            }
        } else {
            subTreeIndex = firstChild(rep->ktree, e.sTree);
            //child function return -1 if sTree no have any child
            //no tiene hijo si en la posición e.sTreee hay un 0.
            if (subTreeIndex != -1) {
                for (int i = 0; i < (K * K); i++) {
                    if (access(rep->ktree,subTreeIndex + i)==1) {
                        //si hay puntos en el subArbol.
                        ElementQueue child;
                        child.sTree = subTreeIndex + i;
                        child.area = getSubArea(e.area, i);
                        child.priority = maxDist(q, child.area);
                        minDis2Child = minDist(q, child.area);
                        if (minDis2Child < pQueue.top().priority) {
                            pQueue.push(child);
                            printPush(child);
                        }
                    }
                }
            }
        }
    }
    //hay que obtener el ID de c.
    if (foundNN) {
        p.setX(c.area.x1);   //en la hoja x1=x2
        p.setY(c.area.y1);   //en la hoja y1=y2
        uint len = 1;
        //TODO: revisar cálculo de rankLeaf, parece que hay un error en c.sTree
        int rankLeaf = rep->ktree->bl->rank1(posInLeaf(rep->ktree, c.sTree));
        uint * resp = rep->labels->getUpToKObjects(rankLeaf, len);
        //está garantizado que hay por lo menos un objeto en la celda por las
        //propiedades del K2-tree.
        oid = resp[0];
        delete[] resp;
        resp = NULL;
    } else {
        oid = -1;
    }

}

//Filtros des vecinos más cercanos.
/*
 * Si la consulta de vecinos más cercanos cae en un instante que no tiene
 * snapshot entonces es necesario realizar la consutalta en el snapshot más
 * cercano.  Esto tiene como consecuencia que no se pueda determinar en el
 * snapshot los vecinos más cercanos, sino que entregar una lista de candidatos
 * posibles, los cuales deben ser evaluados posteriormente, luego de actualizar
 * su posición al instante de la consulta.  Para la realización del filtro
 * tanto en las consultas del vecino más cercano y la de los k vecinos más
 * cercanos se necesita incorporar un parámetro que corresconde a la distancia
 * máxima que un objeto se puede desplazar desde el snapshot al instante de la
 * consulta.  La distancia está dada en número de celdas que es la unidad mínima
 * de distancia en nuestro modelo.
 *
 */

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*
 * oid y p son arreglos no nulos de tamaño kn
 */
void kNNFilter(Snapshot * rep, Point q, uint maxDesp, uint & kn, uint * oid,
                Point * p) {
    priority_queue<CandidateElement, std::vector<CandidateElement>,
            maxHeapComparator> candidates;
    priority_queue<ElementQueue, std::vector<ElementQueue>, minHeapComparator> pQueue;
    ElementQueue e;
    CandidateElement ce;
    e.sTree = 0;
    e.area.x1 = 0;
    e.area.x2 = pow(2, rep->ktree->maxLevel) - 1;
    e.area.y1 = e.area.x1;
    e.area.y2 = e.area.x2;
    e.priority = maxDist(q, e.area);
    long int subTreeIndex;
    double minDis2Child;
    long uint lenCandidates = 0;
    uint numObj = 0;
    uint toAdd = 0;
    uint len = 0;
    uint * resp;
    int rankLeaf;
    pQueue.push(e);
    while (!pQueue.empty()) {
        e = pQueue.top();
        pQueue.pop();
        if (isLeaf(rep->ktree, e.sTree)) {
            //es una hoja de 'ultimo nivel
            //aquí la hoja entera entra a los candidatos o se descarta.
            //por lo tanto la recuperación de los objetos candidatos se deja
            //para el final.

        } else {
            subTreeIndex = firstChild(rep->ktree, e.sTree);
            if (subTreeIndex != -1) {
                //no es una hoja intermedia, es decir una región sin objetos
                for (int i = 0; i < (K * K); i++) {
                    ElementQueue child;
                    child.sTree = subTreeIndex + i;
                    child.area = getSubArea(e.area, i);
                    child.priority = maxDist(q, child.area);
                    minDis2Child = minDist(q, child.area);

                    if (!((candidates.size() >= kn)
                            && (minDis2Child >= candidates.top().priority))) {
                        pQueue.push(child);
                    }
                }
            }
        }
    }
    kn = candidates.size();
    //puede ser que existan menos kn que los originalmente pedidos
    //entonces en la respuesta kn serán los realmente obtenidos.
    //lo siguiente es un overhead que se podría evitar si se retorna la cola
    //directamente.
    uint i = kn - 1;
    while (!candidates.empty()) {
        ce = candidates.top();
        oid[i] = ce.oid;
        p[i] = Point(ce.x, ce.y);
        i--;
        candidates.pop();
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

class compararDistanciasSegmento{
public:
    bool operator()(CHElementQueue& a, CHElementQueue& b){ //El menor valor irá a la cabeza de la cola
        if(a.distancia > b.distancia) return true;
        else if(a.distancia == b.distancia){ //Si tienen la misma distancia, el que tenga menor altura gana
            if(a.h > b.h) return true;
            else return false;
        } else return false;
    }
};

int pendiente(Point a, Point b){
    /*
     * Valores :
     *  1.- Esquina superior izquierda
     *  2.- Esquina superior derecha
     *  3.- Esquina inferior izquierda
     *  4.- Esquina inferior derecha
     * */
    //todo: Revisar viabilidad de punto medio para casos de pendiente 0.
    double pd = 0.0;
    if(b.getX() == a.getX()){ //Verticalmente iguales
        if(b.getY() > a.getY()) return 2;
        else return 1;
    } else if(b.getY() == a.getY()){ //Horizontalmente iguales
        if(b.getX() > a.getX()) return 1;
        else return 3;
    }
    pd = (double)(b.getY() - a.getY())/(b.getX() - a.getX());
    //cout << "PENDIENTE: " << pd << endl;
    if(pd >= 0){
        if(b.getX() > a.getX()) return 2;
        else return 3;
    } else {
        if(b.getX() > a.getX()) return 1;
        else return 4;
    }
    return 0;
}

double giroReloj(Point a, Point b, Point c, uint h, int direccion){ //Punto a y b forman segmento, punto c es el punto a medir contra el segmento
    double gr = 0;
    if(h == 0){
        if((a.getX() == c.getX() && a.getY() == c.getY())||(b.getX()==c.getX() && b.getY()==c.getY())) return 1;
        switch (direccion) {
            case 1:
                if(c.getX() < a.getX() || c.getX() > b.getX() || c.getY() > a.getY() || c.getY() < b.getY()) return 2;
                break;
            case 2:
                if(c.getX() < a.getX() || c.getX() > b.getX() || c.getY() < a.getY() || c.getY() > b.getY()) return 3;
                break;
            case 3:
                if(c.getX() > a.getX() || c.getX() < b.getX() || c.getY() > a.getY() || c.getY() < b.getY()) return 4;
                break;
            case 4:
                if(c.getX() > a.getX() || c.getX() < b.getX() || c.getY() < a.getY() || c.getY() > b.getY()) return 5;
                break;
            default:
                break;
        }
    }

    gr = (b.getX()-a.getX())*(c.getY()-a.getY())-(b.getY()-a.getY())*(c.getX()-a.getX());
    return gr;
}

double giroReloj(Point a, Point b, Point c){ //Punto a y b forman segmento, punto c es el punto a medir contra el segmento
    double gr = 0;
    if((a.getX() == c.getX() && a.getY() == c.getY())||(b.getX()==c.getX() && b.getY()==c.getY())) return 1;
    gr = (b.getX()-a.getX())*(c.getY()-a.getY())-(b.getY()-a.getY())*(c.getX()-a.getX());
    return gr;
}

double distanciaSegmento(Point a, Point b, Point c){ //Los mismos puntos que giroReloj
    double ds = 0.0;
    ds = giroReloj(a, b, c)/sqrt(pow((b.getX()-a.getX()),2)+pow((b.getY()-a.getY()),2));
    return ds;
}

Point definePuntoC(Point p, uint direccion, uint h){
    Point p2;
    switch (direccion) {
        case 1:
            p2.setX(p.getX());
            p2.setY(p.getY());
            break;
        case 2:
            p2.setX(p.getX()+pow(K,h)-1);
            p2.setY(p.getY());
            break;
        case 3:
            p2.setX(p.getX());
            p2.setY(p.getY()+pow(K,h)-1);
            break;
        case 4:
            p2.setX(p.getX()+pow(K,h)-1);
            p2.setY(p.getY()+pow(K,h)-1);
            break;
        default:
            p2 = p;
            break;
    }
    return p2;
}

Point definePuntoC(uint x, uint y, uint direccion, uint h){
    Point p2;
    switch (direccion) {
        case 1:
            p2.setX(x);
            p2.setY(y);
            break;
        case 2:
            p2.setX(x+pow(K,h)-1);
            p2.setY(y);
            break;
        case 3:
            p2.setX(x);
            p2.setY(y+pow(K,h)-1);
            break;
        case 4:
            p2.setX(x+pow(K,h)-1);
            p2.setY(y+pow(K,h)-1);
            break;
        default:
            p2.setX(x);
            p2.setY(y);
            break;
    }
    return p2;
}

bool pertinencia(Point a, Point b, Point c){

    return false;
}

uint setXChild(uint x, uint i, uint h){
    uint av = i%K;
    uint nx = x + av*pow(K,h-1);
    return nx;
}

uint setYChild(uint y, uint i, uint h){
    uint av = floor((i%(K*K))/K);
    uint nx = y + av*pow(K,h-1);
    return nx;
}

void vaciarCola(priority_queue<CHElementQueue, std::vector<chElementQeue>, compararDistanciasSegmento> &pQueue){
    CHElementQueue a;
    while(!pQueue.empty()){
        a = pQueue.top();
        pQueue.pop();
    }
}

Point * desintegra(MREP2 * k2tree, int n){
    Point *p = new Point[n];
    stack<CHElementQueue> pila;
    CHElementQueue aux, auxx;
    aux.sTree = 0;
    aux.x = 0;
    aux.y = 0;
    aux.h = k2tree->maxLevel;
    int indice,cont=0;
    pila.push(aux);
    int subTreeIndex = 0;
    while (!pila.empty()) {
        aux = pila.top();
        pila.pop();
        if (isLeaf(k2tree, aux.sTree)) {
            //si elemento extraído de cola de prioridades es una hoja
            p[cont].setX(aux.x);
            p[cont++].setY(aux.y);
        } else {
            subTreeIndex = firstChild(k2tree, aux.sTree);
            //child function return -1 if sTree no have any child
            //no tiene hijo si en la posición actual.sTreee hay un 0.
            if (subTreeIndex != -1) {
                for (int i = 0; i < 4; i++) {
                    if(access(k2tree, subTreeIndex+i)==0) continue;
                    CHElementQueue child;
                    child.sTree = subTreeIndex + i;
                    child.x = setXChild(aux.x, i, aux.h);
                    child.y = setYChild(aux.y, i, aux.h);
                    child.h = aux.h-1;
                    child.distancia = 0.0;
                    pila.push(child);
                }
            }
        }
    }
    p[cont].setX(-1);
    return p;
}

//stack<Point> convexHull(MREP2 * k2tree, Point n, Point s, Point e, Point o, uint maxlevel) {
int convexHull(MREP2 * k2tree, Point n, Point s, Point e, Point o, uint maxlevel) {
    priority_queue<CHElementQueue, std::vector<CHElementQueue>, compararDistanciasSegmento> pQueue;
    cout << "PQUEUE: " << sizeof(pQueue) << endl;
    //crear una pila con los cuatro puntos extremos.
    stack<Point> * extremos = new stack<Point>();
    //FILE * f2 = fopen("registroConvex.txt","w");
    std::stack<Point> ch;
    extremos->push(Point(-1,-1));
    int flagc = 0;
    extremos->push(o);
    extremos->push(s);
    extremos->push(e);
    extremos->push(n);
    extremos->push(o);
    CHElementQueue actual;
    cout << "Tamano CHElement " << sizeof(CHElementQueue) << endl;
    actual.sTree = 0; //Raíz falsa
    actual.x = 0;
    actual.y = 0;
    actual.h = maxlevel;
    actual.distancia = 0.0;
    long int subTreeIndex;
    pQueue.push(actual);
    uint k2 = K*K;
    Point cand;
    cand.setX(-1);
    //EXTREMOS A ANALIZAR
    Point ex1, ex2;
    ex1 = extremos->top();
    extremos->pop();
    ex2 = extremos->top();
    extremos->pop();
    uint direccion = 0;
    Point c;
    while(ex1.getX() == ex2.getX() && ex1.getY() == ex2.getY()){
        //cout << "EXTREMO PUSH: " << ex1.getX() << "," << ex1.getY() << endl;
    	ch.push(ex1);
        ex1 = ex2;
        ex2 = extremos->top();
        extremos->pop();
    }
    //cout << "EX2: " << ex2.getX() << endl;
    uint tamaCola = 40; //sizeof(CHElementQueue);
    uint max = 0;
    uint maxAux = 0;
    //cout << "no se cae aun" << endl;
    while(!extremos->empty() && ex2.getX() >= 0){
        direccion = pendiente(ex1, ex2);
        while (!pQueue.empty()) {
            actual = pQueue.top();
            pQueue.pop();
            tamaCola -= 40;
            if (isLeaf(k2tree, actual.sTree)) {
                //si elemento extraído de cola de prioridades es una hoja
                cand.setX(actual.x);
                cand.setY(actual.y);
                //if(max < pQueue.size()) max = pQueue.size();
                vaciarCola(pQueue);
            } else {
                subTreeIndex = firstChild(k2tree, actual.sTree);
                //child function return -1 if sTree no have any child
                //no tiene hijo si en la posición actual.sTreee hay un 0.
                if (subTreeIndex != -1) {
                    for (int i = 0; i < k2; i++) {
                        if(access(k2tree, subTreeIndex+i)==0) continue;
                        CHElementQueue child;
                        child.sTree = subTreeIndex + i;
                        child.x = setXChild(actual.x, i, actual.h);
                        child.y = setYChild(actual.y, i, actual.h);
                        child.h = actual.h-1;
                        child.distancia = 0.0;
                        c = definePuntoC(child.x, child.y, direccion, child.h);
                        double gR = giroReloj(ex1, ex2, c, child.h, direccion);
                        //cout << "Cuadrante: " << c.getX() << "," << c.getY() << endl;
                        //cout << "GIRO RELOJ: " << gR << endl;
                        if (gR < 0) {
                            child.distancia = distanciaSegmento(ex1, ex2, c);
                            pQueue.push(child);
                            tamaCola += 40;
                            if(maxAux < tamaCola) maxAux = tamaCola;
                            //cout << "empila" << endl;
                        }
                    }
                }
            }
        }
        //cout << "Candidato: " << cand.getX() << "," << cand.getY() << endl;
        if(cand.getX() >= 0){ //SE OBTUVO UN CANDIDATO
            extremos->push(ex2);
            ex2 = cand;
            cand.setX(-1);
        } else {
            ch.push(ex1);
            ex1 = ex2;
            ex2 = extremos->top();
            extremos->pop();
            while(ex1.getX() == ex2.getX() && ex1.getY() == ex2.getY()){
                ch.push(ex1);
                ex1 = ex2;
                ex2 = extremos->top();
                extremos->pop();
            }
        }
        if(max < maxAux) max = maxAux;
        maxAux = 0;
        tamaCola = 40;
        actual.sTree = 0;
        actual.x = 0;
        actual.y = 0;
        actual.h = maxlevel;
        actual.distancia = 0.0;
        pQueue.push(actual);
    }
    //cout << "TAMANO PQUEUE "<< sizeof(pQueue) << endl;
    //cout << "MAXIMA COLA DE PRIORIDADES: " << max*sizeof(CHElementQueue) << endl;
    return max;
    //return ch;
}

//////////////////////////////////////////////////////////////////////////////
size_t sizeMREP2(MREP2 * rep) {
    size_t totalByte = 0;
    if (rep != NULL) {
        totalByte = sizeof(MREP2);
        //fprintf(stdout,"MREP2: %d\n",sizeof(MREP2));
        totalByte += sizeof(long long uint) * (rep->maxLevel + 1);
        //fprintf(stdout,"TOTAL BYTE: %d\n",sizeof(long long uint) * (rep->maxLevel + 1));
        totalByte += rep->bl->getSize();
        //fprintf(stdout,"BL: %d\n",rep->bl->getSize());
        totalByte += rep->bn->getSize();
        //fprintf(stdout,"BN: %d\n",rep->bn->getSize());
        totalByte += rep->bt->getSize();
        //fprintf(stdout,"BT: %d\n",rep->bt->getSize());
    }
    return totalByte;
}

size_t sizeSnapshot(Snapshot * lktrep) {
    size_t tam = 0;
    if (lktrep != NULL) {
        tam += sizeof(Snapshot);
        //fprintf(stdout,"SNAPSHOT: %d\n",sizeof(Snapshot));
        tam += sizeMREP2(lktrep->ktree);
        //fprintf(stdout,"MREP: %d\n",sizeMREP2(lktrep->ktree));
        //tam += lktrep->labels->getSize();
    }
    return tam;

}
