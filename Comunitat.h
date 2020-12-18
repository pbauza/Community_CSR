#pragma once
#include "MatriuSparse.h"
#include <list>
#include <queue>
#include "Heap.h"
#include <math.h>


#define COMU 0
#define VEI_I 1
#define VEI_J 2

class Comunitat
{
public:
	Comunitat(MatriuSparse* pMAdj);
	~Comunitat() ;
	void calculaM2() { m_M2 = m_pMAdj->getNValues(); }
	void calculaK() { m_pMAdj->calculaGrau(m_k); };
	void calculaA();
	double calcDeltaQ(int i, int j) { return (1 / m_M2) - ((m_k[i] * m_k[j]) / (m_M2 * m_M2)); };
	void creaDeltaQHeap();
	void creaIndexComs();
	void InicialitzaDendrograms(){ m_pMAdj->calculaDendrograms(m_vDendrograms); }
	void calculaComunitats(list<Tree<double>*>& listDendrogram);
	void fusiona(int com1, int com2);
	void modificaVei(int com1, int com2, int vei, int cas);
	int getM2() { return m_M2; }
	vector<int> getK() { return m_k; }
	vector<double> getA() { return m_A; }
	vector<map<pair<int, int>, double>> getdeltaQ() { return m_deltaQ; }
	Heap gethTotal() {return m_hTotal;}
	vector<pair<int, int>> getIndexComs() { return m_indexComs; }
	void clear();
	
private:
	//vector de maps per cada fila de deltaQij
	vector<map<pair<int, int>, double>> m_deltaQ;
	//vector que mante llista de comunitats actives amb un parell que indica anterior i seguent activa.
	//la primera te com anterior -1 i la ultima com a seguent la mida del vector
	vector<pair<int, int>> m_indexComs;
	//vector que mante el maxim de deltaQij d'una fila, per no recalcular innecesariament
	vector<pair<int, double>> m_maxDeltaQFil;
	//index que indica quina es la primera comunitat activa
	int m_primComdeltaQ;
	//Vector d'arbres que va creant el dendograma. Inicialment es un vector amb un arbre per cada node del graf.
	//Cada node fulla te el index del node, cada node intern i l'arrel te el deltaQij que l'ha fet unir-se. Si es negatiu vol dir que no hauriem de fer aquesta unio
	vector<Tree<double>*> m_vDendrograms;
	//Vector de graus dels nodes
	vector<int> m_k;
	//Vector de index aij de cada aresta
	vector<double> m_A;
	//Heap dels maxims deltaQ de cada fila, tindra nfiles elements a l'inici i despres tants com comunitats hi hagi actives
	Heap m_hTotal; // (_mAdjacencia.getNFiles());
	//Modularitat Q de la particio en comunitats calculada
	double m_Q;
	int m_M2; //Técnicamente la matriz es simétrica, así que el número total de valores diferentes a 0 ya es m*2
	MatriuSparse* m_pMAdj;

	int m_nComunitats;
	vector<bool> m_visitats;

};