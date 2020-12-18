#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <ostream>
#include <map>

#include "Tree.hpp"
using namespace std;


class MatriuSparse {
public:

	MatriuSparse() :m_filaPtr(), m_valor(), m_columna(), m_nFiles(0), m_nValors(0), m_nColumnes(0) {} //fet
	~MatriuSparse() { m_filaPtr.clear(); m_valor.clear(); m_columna.clear(); m_nColumnes = 0; m_nFiles = 0; m_nValors = 0; } //fet
	MatriuSparse(int files, int columnes); //fet
	MatriuSparse(const MatriuSparse& m); //fet
	MatriuSparse(string nomFitxer); //fet
	MatriuSparse& operator=(const MatriuSparse &m); //fet
	const bool operator==(const MatriuSparse &m) const; //fet
	MatriuSparse operator*(float v); //fet
	vector<float> operator*(vector<float> v); //fet
	MatriuSparse operator/(float v); //fet
	// void mostra();
	friend ostream& operator<<(ostream& out, const MatriuSparse& M); //fet
	void init(int nFiles, int nColumnes); //fet
	int getNFiles() { return m_nFiles; } //fet
	int getNColumnes() { return m_nColumnes; } //fet
	const bool getVal(const int fila, const int columna, float &valor) const; //fet
	void setVal(const int fila, const int columna, float valor); //fet
	int getNValues() { return m_nValors; }; //Al ser matriz simétrica todos los valores serán m*2
	void calculaGrau(vector<int> &graus) const; //fet
	void calculaDendrograms(vector<Tree<double>*>& vDendograms); //QUEDA HACER ESTO
	void creaMaps(vector<map<pair<int, int>, double>>& vMaps) const; //fet?
	void clear(); //fet

private:

	vector<int> m_filaPtr;
	vector<float> m_valor;
	vector<int> m_columna;
	int m_nFiles;
	int m_nValors;
	int m_nColumnes;
};