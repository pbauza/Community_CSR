#include"MatriuSparse.h"

MatriuSparse::MatriuSparse(int files, int columnes) {
	m_nValors = 0;
	if (files > columnes) {
		m_nFiles = files;
		m_nColumnes = files;
	}
	else {
		m_nFiles = columnes;
		m_nColumnes = columnes;
	}
	m_filaPtr.resize(files, 0);
}

MatriuSparse::MatriuSparse(const MatriuSparse& m) {
	m_nValors = m.m_nValors;
	m_nFiles = m.m_nFiles;
	m_nColumnes = m.m_nColumnes;

	m_filaPtr = m.m_filaPtr;
	m_valor = m.m_valor;
	m_columna = m.m_columna;
}

MatriuSparse& MatriuSparse::operator=(const MatriuSparse& m) {

	if (!(m == *this)) { //falta implementar este operador, por ahora lo dejo asi
		m_nFiles = m.m_nFiles;
		m_nValors = m.m_nValors;
		m_nColumnes = m.m_nColumnes;

		m_filaPtr = m.m_filaPtr;
		m_valor = m.m_valor;
		m_columna = m.m_columna;

	}

	return *this;
}

const bool MatriuSparse::operator==(const MatriuSparse & m) const
{
	bool Returner = false;

	if ((m_nFiles == m.m_nFiles) && (m_nColumnes == m.m_nColumnes) && (m_nValors == m.m_nValors) && (m_filaPtr == m.m_filaPtr) && (m_columna == m.m_columna) && (m_valor == m.m_valor))
		Returner = true;

	return Returner;
}

MatriuSparse MatriuSparse::operator*(float v) {
	MatriuSparse m(*this);

	for (int i = 0; i < m_valor.size(); i++) {
		m.m_valor[i] *= v;
	}

	return m;
}

vector<float> MatriuSparse::operator*(vector<float> v) {
	try {
		if (v.size() == m_nFiles)
		{
			vector<float> VAux;

			VAux.resize(v.size(), 0);

			int counter;
			float valor;
			bool Stop = false;

			for (int i = 0; i < v.size(); i++) {
				counter = m_filaPtr[i];

				while ((counter < m_filaPtr[i + 1]) && (!Stop))
				{
					VAux[i] += (v[i] * m_valor[counter]);

					counter++;
				}
			}

			return VAux;
		}
		else
			throw "NUMERO DE FILES INCORRECTE";
	}
	catch (const char* msg) {
		cout << msg << endl;
	}
}

MatriuSparse MatriuSparse::operator/(float v) {
	try {
		if (v != 0)
		{
			MatriuSparse m(*this);

			for (int i = 0; i < m_nValors; i++) {
				m.m_valor[i] /= v;
			}

			return m;
		}
		else
			throw "NO ES POT DIVIDIR ENTRE 0";
	}
	catch (const char* msg) {
		cout << msg << endl;
	}
}

void MatriuSparse::init(int files, int columnes) {
	m_nValors = 0;
	m_nFiles = files;
	m_nColumnes = columnes;
	m_valor.clear();
	m_columna.clear();
	m_filaPtr.clear();
	m_filaPtr.resize(files);
	for (int i = 0; i < m_nFiles; i++) {
		m_filaPtr[i] = 0;
	}
}

MatriuSparse::MatriuSparse(string nomFitxer) {
	std::ifstream file;
	m_nValors = 0;

	file.open(nomFitxer);
	if (file.is_open()) {
		int fila, columna, fila_aux, fila_aux_1, fila_max, columna_max;
		file >> fila;
		fila_aux = fila;
		fila_max = fila;
		file >> columna;
		columna_max = columna;

		m_valor.resize(m_nValors + 1);
		m_valor[m_nValors] = 1;

		m_columna.resize(m_nValors + 1);
		m_columna[m_nValors] = columna;

		m_filaPtr.resize(fila + 1);
		m_filaPtr[fila] = m_nValors;

		m_nValors++;

		while (!file.eof()) {
			int fila_anterior = fila;
			int columna_anterior = columna;
			file >> fila;
			file >> columna;


			if ((fila != fila_anterior) || (columna != columna_anterior)) {
				if (fila > fila_max) {
					fila_max = fila;
				}
				if (columna > columna_max) {
					columna_max = columna;
				}

				m_valor.resize(m_nValors + 1);
				m_valor[m_nValors] = 1;

				m_columna.resize(m_nValors + 1);
				m_columna[m_nValors] = columna;

				if (fila_aux != fila) {

					for (int filaCtr = fila_aux + 1; filaCtr <= fila; filaCtr++) {
						m_filaPtr.resize(filaCtr + 1);
						m_filaPtr[filaCtr] = m_nValors;
					}

					m_filaPtr[fila] = m_nValors;
					fila_aux_1 = fila_aux;
				}
				fila_aux = fila;
				m_nValors++;
			}
		}

		file.close();

		if (fila_aux_1 != fila_aux) { //caso ultimo valor
			fila++;
			m_filaPtr.resize(fila + 1);
			m_filaPtr[fila] = m_nValors;
		}
		if (fila_max > columna_max) {
			m_nFiles = fila_max + 1;
			m_nColumnes = fila_max + 1;
		}
		else {
			m_nFiles = columna_max + 1;
			m_nColumnes = columna_max + 1;
			m_filaPtr.resize(m_nFiles + 1);
			for (int i = fila_max + 1; i < m_filaPtr.size(); i++) { //meter valores (0) en las filas
				m_filaPtr[i] = m_nValors;
			}
		}
		//las matrices seran siempre cuadradas
	}
}

ostream& operator<<(ostream& out, const MatriuSparse& M) {
	out << "MATRIU DE FILES: " << M.m_nFiles << " : " << "COLUMNES: " << M.m_nColumnes << endl;
	float valor;
	for (int i = 0; i < M.m_filaPtr.size() - 1; i++) {
		if (M.m_filaPtr[i + 1] - M.m_filaPtr[i] > 0) {
			out << "VALORS FILA:" << i << "(COL:VALOR)" << endl;
			int counter = M.m_filaPtr[i];

			if (i < M.m_filaPtr.size()) {
				while ((counter < M.m_filaPtr[i + 1]))
				{
					M.getVal(i, M.m_columna[counter], valor);

					if (valor != 0)
					{
						out << "(" << M.m_columna[counter] << " : " << M.m_valor[counter] << ") ";

						counter++;
					}
				}
				out << endl;
			}
		}

	}

	out << "MATRIUS" << endl << "VALORS" << endl << "(";
	for (int i = 0; i < M.m_nValors; i++)
		out << M.m_valor[i] << "  ";
	out << ")" << endl;
	out << "COLS" << endl << "(";
	for (int i = 0; i < M.m_nValors; i++)
		out << M.m_columna[i] << "  ";
	out << ")" << endl;

	out << "INIFILA" << endl;
	out << "(";
	for (int i = 0; i < M.m_nFiles; i++)
		if (M.m_filaPtr[i] != M.m_filaPtr[i + 1]) {
			out << "[ " << i << " : " << M.m_filaPtr[i] << " ] ";
		}

	out << " [Num Elems:" << M.m_nValors << "] )" << endl;

	return out;
}

const bool MatriuSparse::getVal(const int fila, const int columna, float& valor) const {
	bool Stop = false;
	bool valid = false;

	if ((fila >= 0) && (columna >= 0) && (fila <= m_nFiles - 1) && (columna <= m_nColumnes - 1))
	{
		int counter;
		valid = true;
		valor = 0;

		counter = m_filaPtr[fila];

		while ((counter < m_filaPtr[fila + 1]) && (!Stop))
		{
			if (m_columna[counter] == columna)
			{
				valor = m_valor[counter];
				Stop = true;
			}
			else
				counter++;
		}
	}

	return valid;
}

void MatriuSparse::setVal(const int fila, const int columna, float valor)
{
	int fila_aux = m_nFiles;
	if ((valor != 0) && (fila >= 0) && (columna >= 0))
	{

		float valor_2;
		bool valid = getVal(fila, columna, valor_2);

		if ((fila > m_nFiles) || (columna > m_nColumnes))
		{
			if (fila > columna)
			{
				m_nFiles = fila + 1; //+1 porque hay que tener en cuenta la posicion 0
				m_nColumnes = fila + 1;
				m_filaPtr.resize(fila + 1);
			}
			else
			{
				m_nFiles = columna + 1;
				m_nColumnes = columna + 1;
				m_filaPtr.resize(columna + 1);
			}
		}


		if (!valid || valor_2 == 0) {
			int counter;
			if (fila >= fila_aux) {
				counter = fila_aux;
				m_filaPtr.resize(fila + 2);
				for (counter; counter <= fila + 1; counter++)
				{
					m_filaPtr[counter] = m_nValors;
				}
			}
			else {
				counter = fila + 1;
				m_filaPtr.resize(m_nFiles + 1);
				for (counter; counter <= m_nFiles; counter++)
				{
					m_filaPtr[counter]++;
				}
			}

			m_nValors++;

			//m_valor.resize(m_nValors);
			//m_columna.resize(m_nValors);
			//m_filaPtr.resize(m_filaPtr.size() + 1);

			m_filaPtr[m_filaPtr.size() - 1] = m_nValors;


			vector<float>::iterator it_v = m_valor.begin();
			vector<int>::iterator it_c = m_columna.begin();

			for (int i = 0; i < m_filaPtr[fila]; i++) {
				it_v++;
				it_c++;
			}

			m_valor.insert(it_v, valor);
			m_columna.insert(it_c, columna);
			int aux_c;
			float aux_v;

			vector<int> columna_auxiliar_v;
			columna_auxiliar_v = m_columna;
			int pos = 0;

			for (int i = m_filaPtr[fila]; i < m_filaPtr[fila + 1]; i++) {
				for (int x = i + 1; x < m_filaPtr[fila + 1]; x++) { //ojo
					if (m_columna[i] > m_columna[x]) {
						aux_c = m_columna[i];
						aux_v = m_valor[i];
						m_columna[i] = m_columna[x];
						m_valor[i] = m_valor[x];
						m_columna[x] = aux_c;
						m_valor[x] = aux_v;
					}
				}
			}
		}
		else {
			int i;
			bool Stop = false;

			i = m_filaPtr[fila];

			while ((i < m_filaPtr[fila + 1]) && (!Stop))
			{
				if (m_columna[i] == columna)
				{
					m_valor[i] = valor;
					Stop = true;
				}
				else
					i++;
			}
		}
	}
}


//.-----------------PRACTICA 2------------------.\\

void MatriuSparse::calculaGrau(vector<int> &graus) const
{
	graus.resize(m_filaPtr.size() - 1);
	
	for (int i = 0; i < graus.size(); i++)
		graus[i] = m_filaPtr[i + 1] - m_filaPtr[i];
}

void MatriuSparse::creaMaps(vector<map<pair<int, int>, double>>& vMaps) const {
	vector<int> vGrau;
	map<pair<int, int>, double> mapa;

	this->calculaGrau(vGrau);
	
	int counter = 0;
	int j = 0;
	vMaps.resize(m_filaPtr.size() - 1);

	for (int i = 0; i < vMaps.size(); i++) {
		for (j = counter; j < counter + vGrau[i]; j++) {
			pair<int, int> parella;
			parella.first = i;
			parella.second = m_columna[j];
			
			mapa.insert(pair<pair<int, int>, double>(parella, 0));
		}
		counter = j;
		vMaps[i] = mapa;
		mapa.clear();
	}
}

void MatriuSparse::calculaDendrograms(vector<Tree<double>*>& vDendograms) { //al fin, no tocar que esta bien

	for (double i = 0; i < m_filaPtr.size() - 1; i++) {
		vDendograms.resize(i + 1);
		vDendograms[i] = new Tree<double>(i);
	}
}

void MatriuSparse::clear() {
	m_columna.clear();
	m_filaPtr.clear();
	m_valor.clear();
	m_nColumnes = 0;
	m_nFiles = 0;
	m_nValors = 0;
}