#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>
#include <utility>
#include <iomanip>
#include "/mnt/c/Users/dravv/Projects/Qt/BodeDiagram/exprtk/exprtk.hpp" 

class MagnitudeAndPhase {
    private:
        std::string _userFunction;
        int _freqMin;
        int _freqMax;
	double _s_real;

    public:
        MagnitudeAndPhase(std::string userFunction, int freqMin, int freqMax, double s_real );
        std::pair<std::vector<double>, std::vector<double>> frequencies(); 
        std::string translatefunction(double magnitude);
	double evaluatetranslatedfunction(std::string translated);
        double calculateMagnitude(double w);
        double calculatePhase(double w);
};

// Constructor
MagnitudeAndPhase::MagnitudeAndPhase(std::string userFunction, int freqMin, int freqMax, double s_real) 
    : _userFunction(userFunction), _freqMin(freqMin), _freqMax(freqMax), _s_real(s_real){}

// Definición de la función frequencias
std::pair<std::vector<double>, std::vector<double>> MagnitudeAndPhase::frequencies() {
    std::vector<double> freqVector;
    std::vector<double> angularFreqVector;

    for (double i = _freqMin; i <= _freqMax; i += 1.0) {
        freqVector.push_back(i);
        angularFreqVector.push_back(2 * M_PI * i); // Cálculo de la frecuencia angular
    }

    return std::make_pair(angularFreqVector, freqVector); // Retorna el par de vectores
}
//Definicion de la funcion de transferir valores de s
std::string MagnitudeAndPhase::translatefunction(double magnitude){
	std::string translatedfunction = _userFunction;
	size_t pos = 0;

	std::string magnitudeStr = std::to_string(magnitude);
	while ((pos = translatedfunction.find("s", pos)) != std::string::npos) {
		translatedfunction.replace(pos, 1, magnitudeStr);
		pos += magnitudeStr.length();
	}
	return translatedfunction;
}

// Calcula la magnitud |s| = |j * w|
double MagnitudeAndPhase::calculateMagnitude(double w) {
    std::complex<double> s(_s_real, w);  // s = j * w
    return std::abs(s);  // Retorna la magnitud |s|
}
double MagnitudeAndPhase::evaluatetranslatedfunction(std::string translated){
	exprtk::expression<double> expression;
	exprtk::parser<double> parser;

	if (!parser.compile(translated, expression)) {
		std::cerr << "Error al compilar la funcion: "<<parser.error() << std::endl;
		return 0.0; 
	}
	double result = expression.value();
	double magnitudedB = 20 *log10(result);
	return magnitudedB; 

}
double MagnitudeAndPhase::calculatePhase(double w){
	std::complex<double> s(_s_real, w);
        double magnitude = std::abs(s);
        double phase = std::arg(s);

        std::cout << "s: " << s << " -> Magnitud: " << magnitude << ", Fase en (radianes): " << phase << std::endl;

	return std::arg(s)*(180.0 /M_PI);//retorna la fase en angulo 
}
//tancan(s+1)


int main() {
    double _freqMin;
    double _freqMax;
    std::string _userFunction;
    double _s_real;

    std::cout << "escribe la funcion: " << std::endl;
    std::getline(std::cin, _userFunction);
    std::cout << "escribe la freqMin: " << std::endl;
    std::cin >> _freqMin;
    std::cout << "escribe la freqMax: " << std::endl;
    std::cin >> _freqMax;
    std::cout << "escribe el valor de sigma (o) si es diferente a 0.0: "<< std::endl;
    std::cin >> _s_real;

    // Crear objeto de la clase MagnitudeAndPhase
    MagnitudeAndPhase mapObject(_userFunction, _freqMin, _freqMax, _s_real);

    //llama a la funcion frequencies();
    auto [angularFreqVector, freqVector] = mapObject.frequencies();

    // Imprimir frecuencias angulares
    std::cout << "frequencias angulares Totales (w): ";
    for (const auto& w : angularFreqVector) {
        std::cout << std::fixed << std::setprecision(4) << w << " "; // Agregar un espacio entre los valores
    }
    std::cout << std::endl;

    // Imprimir frecuencias en Hz
    std::cout << "frequencias Totales (Hz): ";
    for (const auto& freq : freqVector) {
        std::cout << std::fixed << std::setprecision(4) << freq << " "; // Agregar un espacio entre los valores
    }
    std::cout << std::endl;


    //resultado de funcion traducida translatedfunction() y
    // Calcular e imprimir la magnitud de s para cada w
    std::cout << "Magnitud y Fase: ";
    for (const auto& w : angularFreqVector) {
        double magnitude = mapObject.calculateMagnitude(w);
	double phase = mapObject.calculatePhase(w);

	std::cout << " --> Magnitud: "<< std::fixed << std::setprecision(4) << magnitude << std::endl;
	std::cout <<"--> Fase Total en (grados): "<< std::fixed << std::setprecision(4) << phase << std::endl;
        

    //llamar a la funcion translatedfunction(con magnitude incluida)
    std::cout <<"Formula (despues de traducir): ";
    std::string translated = mapObject.translatefunction(magnitude);
    std::cout <<"--> " << translated << " ";


    //llama a evaluatetranslatedfunction(translated) para evaluar
    double magnitudedB = mapObject.evaluatetranslatedfunction (translated);
    std::cout << "Magnitud Total en (dB): "<< std::fixed<< std::setprecision(4) << magnitudedB << std::endl;
    }
    std::cout << std::endl;




    return 0;
}

