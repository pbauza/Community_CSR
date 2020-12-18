#include "Comunitat.h"

Comunitat::Comunitat(MatriuSparse* pMAdj)
{
	m_primComdeltaQ = -1;
	m_Q = 0;
	m_pMAdj = pMAdj;
}

Comunitat::~Comunitat()
{
	m_pMAdj = nullptr;	
}

void Comunitat::clear()
{
	m_pMAdj = nullptr;
	m_deltaQ.clear();
	m_indexComs.clear();
	m_maxDeltaQFil.clear();
	m_primComdeltaQ = -1;
	m_vDendrograms.clear();
	m_k.clear();
	m_A.clear();
	m_hTotal.clear();
	m_Q = 0;
	m_M2 = 0;
}

void Comunitat::calculaA()
{
	m_A.resize(m_k.size(), 2);

	for (int i = 0; i < m_A.size(); i++)
		m_A[i] = double(m_k[i]) / double(m_M2);
}

void Comunitat::creaIndexComs()
{
	m_primComdeltaQ = 0;
	m_nComunitats = m_k.size();
	m_visitats.resize(m_k.size(), false);

	m_indexComs.resize(m_k.size());
	for (int i = 0; i < m_indexComs.size(); i++) {
		m_indexComs[i].second = i - 1; //deberia estar al reves
		m_indexComs[i].first = i + 1;
	}
}

void Comunitat::creaDeltaQHeap() //este lo dejamos?? lo del heap al final?
{
	m_pMAdj->creaMaps(m_deltaQ);
	m_hTotal.resize(m_k.size());
	pair<int, double> max_p;
	map<pair<int, int>, double>::iterator it;
	double max;
	m_maxDeltaQFil.resize(m_k.size());
	
	for (int i = 0; i < m_deltaQ.size(); i++) {
		it = m_deltaQ[i].begin();
		max = -2;
		int nodeAdj;
		for (int j = 0; j < m_deltaQ[i].size(); j++) {
			nodeAdj = it->first.second;
			it->second = (1 / double(m_M2)) - ((double(m_k[i]) *double(m_k[nodeAdj])) / (double(m_M2) * double(m_M2)));
			if (max < it->second) {
				max = it->second;
				//max = roundf(max * pow(10,7)) / pow(10, 7);
				max_p.first = it->first.second;
				max_p.second = max;
			}

			it++;
		}

		m_maxDeltaQFil[i] = max_p;
		ElemHeap aInserir;
		pair<int, int> pos;
		pos.first = i;
		//pos.second = max_p.first;

		if (m_deltaQ[i].empty()) {
			pos.second = -1;
			//max = 0;
		}
		else {
			pos.second = max_p.first;
		}

		aInserir = ElemHeap(max, pos);
		m_hTotal.insert(aInserir);
	}
}

void Comunitat::modificaVei(int com1, int com2, int vei, int cas)
{
	if (cas == COMU) {
		m_deltaQ[vei].at({ vei, com2 }) += m_deltaQ[vei].at({ vei, com1 });
		m_deltaQ[com2].at({ com2, vei }) += m_deltaQ[com1].at({ com1, vei });
		m_deltaQ[vei].erase({ vei, com1 });
		m_maxDeltaQFil[vei].second = -2; //el posem a -2 perque despres es modifiqui si o si
		bool modificat = false;

		map<pair<int, int>, double>::iterator itvei = m_deltaQ[vei].begin();
		map<pair<int, int>, double>::iterator itaux;
		while (itvei != m_deltaQ[vei].end()) {
			if (itvei->second > m_maxDeltaQFil[vei].second) { //si el deltaq es major, modifiquem el max
				modificat = true;
				m_maxDeltaQFil[vei].second = itvei->second; //modifiquem
				m_maxDeltaQFil[vei].first = itvei->first.second; //modifiquem relacio
				itaux = itvei;
			}
			itvei++;
		}
		if (modificat)
			m_hTotal.modifElem(ElemHeap(itaux->second, { vei, itaux->first.second }));

		//MODIFICAR MAXDELTA I HEAP DE LA J
		if (m_maxDeltaQFil[com2].second < m_deltaQ[com2].at({ com2, vei })) { //si el maxim es menor que el deltaq amb el vei
			m_maxDeltaQFil[com2].second = m_deltaQ[com2].at({ com2, vei }); //modifiquem el max
			m_maxDeltaQFil[com2].first =vei; //modifiquem relacio
			m_hTotal.modifElem(ElemHeap(m_deltaQ[com2].at({ com2, vei }), { com2, vei }));
		}
	}

	if (cas == VEI_I) {
		double valor = m_deltaQ[vei].at({ vei, com1 }) - 2 * double(m_A[com2]) * double(m_A[vei]);
		double valor1 = m_deltaQ[com1].at({ com1, vei }) - 2 * double(m_A[com2]) * double(m_A[vei]);
		m_deltaQ[vei].insert(pair<pair<int, int>, double> { {vei, com2}, valor});
		m_deltaQ[vei].erase({ vei, com1 });
		//m_maxDeltaQFil[vei].second = -2; //el posem a -2 perque despres es modifiqui si o si
		bool modificat = false;

		map<pair<int, int>, double>::iterator itvei = m_deltaQ[vei].begin();
		map<pair<int, int>, double>::iterator itaux;
		while (itvei != m_deltaQ[vei].end()) {
			if (itvei->second > m_maxDeltaQFil[vei].second) { //si el deltaq es major, modifiquem el max
				modificat = true;
				m_maxDeltaQFil[vei].second = itvei->second; //modifiquem
				m_maxDeltaQFil[vei].first = itvei->first.second; //modifiquem relacio
				itaux = itvei;
			}
			itvei++;
		}
		if (modificat)
			m_hTotal.modifElem(ElemHeap(itaux->second, { vei, itaux->first.second }));

		m_deltaQ[com2].insert(pair<pair<int, int>, double> { {com2, vei}, valor1});

		//MODIFICAR MAXDELTA I HEAP DE LA J
		if (m_maxDeltaQFil[com2].second < m_deltaQ[com2].at({ com2, vei })) { //si el maxim es menor que el deltaq amb el vei
			m_maxDeltaQFil[com2].second = m_deltaQ[com2].at({ com2, vei }); //modifiquem el max
			m_maxDeltaQFil[com2].first = vei; //modifiquem relacio
			m_hTotal.modifElem(ElemHeap(m_deltaQ[com2].at({ com2, vei }), { com2, vei }));
		}
	}

	if (cas == VEI_J) {
		m_deltaQ[vei].at({ vei, com2 }) -= 2 * double(m_A[com1]) * double(m_A[vei]);
		m_deltaQ[com2].at({ com2, vei }) -= 2 * double(m_A[com1]) * double(m_A[vei]);
		bool modificat = false;

		map<pair<int, int>, double>::iterator itvei = m_deltaQ[vei].begin();
		map<pair<int, int>, double>::iterator itaux;
		while (itvei != m_deltaQ[vei].end()) {
			if (itvei->second > m_maxDeltaQFil[vei].second) { //si el deltaq es major, modifiquem el max
				modificat = true;
				m_maxDeltaQFil[vei].second = itvei->second; //modifiquem
				m_maxDeltaQFil[vei].first = itvei->first.second; //modifiquem relacio
				itaux = itvei;
			}
			itvei++;
		}
		if(modificat)
			m_hTotal.modifElem(ElemHeap(itaux->second, { vei, itaux->first.second }));
	}


}

void Comunitat::fusiona(int com1, int com2)
{
	//map<pair<int, int>, double> pendents;
	//vector<int> tipus_vei;
	//tipus_vei.resize(m_k.size(), -1);
	map<pair<int, int>, double>::iterator itcom2 = m_deltaQ[com2].begin();
	
	int vei;
	bool Stop;

	//m_hTotal.modifElem(ElemHeap(-2, { com1, com2 }));

	while (itcom2 != m_deltaQ[com2].end()) {
		Stop = false;
		map<pair<int, int>, double>::iterator itcom1 = m_deltaQ[com1].begin();
		vei = itcom2->first.second;
		while (itcom1 != m_deltaQ[com1].end() && !Stop) {
			//pendents.insert({ {com1, itcom2->first.second}, itcom2->second });

			//MODIFIQUEM ELS VEINS EN COMU

			if (vei == itcom1->first.second) { //veins en comu
				//pendents.erase({ com1, itcom2->first.second });
				//tipus_vei[vei] = COMU;
				modificaVei(com1, com2, vei, COMU);
				Stop = true;
			}
			itcom1++;

			//-----------------------------
		}

		if (!Stop) { //vei nomes de la j

			//MODIFIQUEM ELS VEINS NOMES DE LA J

			if (!Stop && vei != com2) { //vei nomes de la j
				//pendents.erase({ com1, itcom2->first.second });
				//tipus_vei[vei] = VEI_J;
				modificaVei(com1, com2, vei, VEI_J);
			}

			//------------------------------------
		}
		itcom2++;
	}

	//---------------------------------------------------------------



	//MODIFIQUEM VEINS NOMES DE LA COMUNITAT 1 (LA I)

	//map<pair<int, int>, double>::iterator it = pendents.begin();
	map<pair<int, int>, double>::iterator it_i = m_deltaQ[com1].begin();
	


	while (it_i != m_deltaQ[com1].end()) {
		bool trobat = false;
		map<pair<int, int>, double>::iterator it_j = m_deltaQ[com2].begin();
		while (it_j != m_deltaQ[com2].end() && !trobat) {
			if (it_i->first.second == it_j->first.second) {
				trobat = true;
			}
			it_j++;
		}
		if (!trobat && it_i->first.second != com2)
			modificaVei(com1, com2, it_i->first.second, VEI_I);
		it_i++;
	}

	//---------------------------------------------------


	// MODIFIQUEM EL MAXDELTA I EL HEAP DE LA COMUNITAT 2 SI EL SEU MAXIM RESIDEIX EN LA RELACIO COM2 - COM1

	map<pair<int, int>, double>::iterator auxx = m_deltaQ[com2].find({ com2, com1 });

	if (auxx != m_deltaQ[com2].end()) {
		if (m_maxDeltaQFil[com2].first == auxx->first.second)
		{
			m_deltaQ[com2].erase(auxx); //ELIMINEM LA RELACIO COM2 - COM1
			double max = -2;
			int i = com2, j = -1;
			map<pair<int, int>, double>::iterator its;
			for (its = m_deltaQ[com2].begin(); its != m_deltaQ[com2].end(); its++)
			{
				if (its->second > max)
				{
					max = its->second;
					i = its->first.first;
					j = its->first.second;
				}
			}
			m_maxDeltaQFil[com2].first = j;
			m_maxDeltaQFil[com2].second = max;
			ElemHeap nou(max, { i, j });
			m_hTotal.modifElem(nou);
		}
		else
			m_deltaQ[com2].erase({ com2, com1 }); //ELIMINEM LA RELACIO COM2 - COM1
	}	

	m_deltaQ[com1].clear(); //eliminem tots els veins de la i

	//-----------------------------------------------------------------------------------



	//ACTUALITZEM LES COMUNITATS ACTIVES/ELIMINADES
	if (m_indexComs[com1].first < m_indexComs.size())
		m_indexComs[m_indexComs[com1].first].second = m_indexComs[com1].second;
	if(m_indexComs[com1].second > -1)
		m_indexComs[m_indexComs[com1].second].first = m_indexComs[com1].first;

	if (com1 == 0) {
		m_primComdeltaQ = m_indexComs[m_primComdeltaQ].first; //si esborrem la primera, augmentem m_primcomdelta
	}
	//---------------------------------------------



	// Creo que esta parte de los arboles esta bien asi

	Tree<double>* arbreIns;
	arbreIns = new Tree<double>(m_Q);

	arbreIns->setLeft(m_vDendrograms[com2]);
	arbreIns->setRight(m_vDendrograms[com1]);

	m_vDendrograms[com2] = arbreIns;
	m_vDendrograms[com1] = NULL;

	m_A[com2] += m_A[com1];
	m_nComunitats--;
	m_visitats[com1] = true;
}

void Comunitat::calculaComunitats(list<Tree<double>*>& listDendrogram)
{
	calculaM2();
	calculaK();
	calculaA();
	creaDeltaQHeap();
	creaIndexComs();
	InicialitzaDendrograms();
	m_Q = 0;

	while (m_nComunitats > 1 && m_hTotal.max().getVal() > 0) {
		if (m_hTotal.max().getVal() + m_Q > m_Q) {
			int com1 = m_hTotal.max().getPos().first;
			int com2 = m_hTotal.max().getPos().second;
			m_hTotal.delMax();
			if (com1 != com2) {
				m_Q += m_hTotal.max().getVal();
				fusiona(com1, com2);
			}
		}else
			m_hTotal.delMax();
	}

	for (int i = m_primComdeltaQ; i < m_vDendrograms.size(); i = m_indexComs[i].first) {
			listDendrogram.push_back(m_vDendrograms[i]);
	}

}

