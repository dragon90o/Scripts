#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>
#include <regex>
#include <utility>
#include <iomanip>
#include "/mnt/c/Users/dravv/Projects/Qt/BodeDiagram/exprtk/exprtk.hpp" 
#include "eigen-3.4.0/eigen-3.4.0/Eigen/Dense"

class MagnitudeAndPhase {
    private:
        std::string _numerator;
	std::string _denominator;
        int _freqMin;
        int _freqMax;
	double _s_real;
	std::vector<std::complex<double>> poles;
	std::vector<std::complex<double>> zeros;

    public:
        MagnitudeAndPhase(std::string numerator, std::string denominator, int freqMin, int freqMax, double s_real );
	std::vector<double> parseCoefficients(const std::string& polynomial);
	std::vector<std::complex<double>> findRoots(const std::vector<double>& coefficients);
        std::pair<std::vector<double>, std::vector<double>> frequencies(); 
        std::string translatefunction(double magnitude, bool isNumerator);
	double evaluatetranslatedfunction(double magnitude);
        double calculateMagnitude(double w);
        double calculatePhase(double w);
	void processTransferFunction();

};

// Constructor fuera de la clase.
MagnitudeAndPhase::MagnitudeAndPhase(std::string numerator, std::string denominator, int freqMin, int freqMax, double s_real) 
    : _numerator(numerator), _denominator(denominator), _freqMin(freqMin), _freqMax(freqMax), _s_real(s_real){}

//Metodo: convierte la cadena de un polinomio en una lista de coeficientes numericos. 
//-- polynomial es la estancia para pasar  numerador y denominador 
std::vector<double> MagnitudeAndPhase::parseCoefficients(const std::string& polynomial) {
    std::vector<double> coefficients;
    std::string modifiedPoly = polynomial;

    
    // Usar una expresión regular para separar la expresión y manejar correctamente los operadores.
    std::regex termRegex("([+-]?\\d*\\.\\d+|[+-]?\\d+)(s\\^\\d+|s)?"); // Regex para manejar coeficientes con s^n o s
    auto termsBegin = std::sregex_iterator(modifiedPoly.begin(), modifiedPoly.end(), termRegex);
    auto termsEnd = std::sregex_iterator();

    // Iterar sobre los términos y extraer los coeficientes
    for (std::sregex_iterator i = termsBegin; i != termsEnd; ++i) {
        std::smatch match = *i;
        std::string coefStr = match[1].str();  // El coeficiente (con signo, si existe)
        std::string sTerm = match[2].str();  // El término 's' o 's^n'

        // Convertir el coeficiente a número
        try {
            double coefficient = std::stod(coefStr);
            coefficients.push_back(coefficient);
        } catch (...) {
            std::cerr << "Error al convertir el término: " << coefStr << std::endl;
        }
    }


    return coefficients;
}
    
//Metodo: calcula las raices (polos y ceros) de un polinomio dado sus coeficientes.
//--
std::vector<std::complex<double>> MagnitudeAndPhase::findRoots(const std::vector<double>& coefficients) {
    int degree = coefficients.size() - 1;
    Eigen::MatrixXd companion(degree, degree);
    companion.setZero();

    for (int i = 0; i < degree; ++i) {
        companion(i, degree - 1) = -coefficients[i] / coefficients.back();
        if (i < degree - 1) {
            companion(i + 1, i) = 1.0;
        }
    }

    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
    Eigen::VectorXcd roots = solver.eigenvalues();

    std::vector<std::complex<double>> result(roots.size());
    for (int i = 0; i < roots.size(); ++i) {
        result[i] = roots[i];
    }
    return result;
}
//Metodo: Procesa la función de transferencia: calcula ceros, polos e imprime los resultados.
//--
void MagnitudeAndPhase::processTransferFunction() {
    // Convertir cadenas a coeficientes
    auto numCoefficients = parseCoefficients(_numerator);
    auto denCoefficients = parseCoefficients(_denominator);

    // Calcular ceros y polos
    zeros = findRoots(numCoefficients);
    poles = findRoots(denCoefficients);

    // Imprimir resultados
    std::cout << "Ceros (Raíces del Numerador): ";
    for (const auto& zero : zeros) {
        std::cout << zero << " ";
    }
    std::cout << std::endl;

    std::cout << "Polos (Raíces del Denominador): ";
    for (const auto& pole : poles) {
        std::cout << pole << " ";
    }
    std::cout << std::endl;
}

//Metodo: genera dos vectores; frecuencias en hz y frecuencias angulares (w).
//--
std::pair<std::vector<double>, std::vector<double>> MagnitudeAndPhase::frequencies() {
    std::vector<double> freqVector;
    std::vector<double> angularFreqVector;

    for (double i = _freqMin; i <= _freqMax; i += 1.0) {
        freqVector.push_back(i);
        angularFreqVector.push_back(2 * M_PI * i); // Cálculo de la frecuencia angular
    }

    return std::make_pair(angularFreqVector, freqVector); // Retorna el par de vectores
}
//Metodo: Remplaza las variables `s` en el numerador o denominador con el valor de magnitud.
//--
std::string MagnitudeAndPhase::translatefunction(double magnitude, bool isNumerator) {
    std::string translatedfunction; 
    if (isNumerator) {
        translatedfunction = _numerator;  // Si estamos procesando el numerador
    } else {
        translatedfunction = _denominator;  // Si estamos procesando el denominador
    }

     // Reemplazar todos los 's' por el valor de la magnitud
    std::string magnitudeStr = std::to_string(magnitude);
    size_t pos = 0;
    while ((pos = translatedfunction.find("s", pos)) != std::string::npos) {
        translatedfunction.replace(pos, 1, "(" + magnitudeStr + ")");
        pos += magnitudeStr.length();
    }

    return translatedfunction;
}

//Metodo: Calcula la magnitud s ," |s| = |j * w| " para una frecuencia angular dada.
//--
double MagnitudeAndPhase::calculateMagnitude(double w) {
    
    std::complex<double> s(_s_real, w);  // s = j * w
    return std::abs(s);  // Retorna la magnitud |s|
}

//Metodo: Evalúa la función de transferencia con la magnitud dada (numerador / denominador).
//--
double MagnitudeAndPhase::evaluatetranslatedfunction(double magnitude){

        std::string translatedNumerator = translatefunction(magnitude, true);
        std::string translatedDenominator = translatefunction(magnitude, false);

	//evalua numerador
	exprtk::expression<double> expressionNumerator;
	exprtk::parser<double> parserNumerator;
	if (!parserNumerator.compile(translatedNumerator, expressionNumerator)) {
		std::cerr << "Error al compilar la funcion: "<<parserNumerator.error() << std::endl;
		return 0.0; 
	}
	double resultNumerator = expressionNumerator.value();

	//evalua Denominador
	exprtk::expression<double> expressionDenominator;
	exprtk::parser<double> parserDenominator;
	if (!parserDenominator.compile(translatedDenominator, expressionDenominator)) {
        std::cerr << "Error al compilar el denominador: " << parserDenominator.error() << std::endl;
        return 0.0;
    }
	double resultDenominator = expressionDenominator.value();

	double transferFunction= resultNumerator / resultDenominator;

	double magnitudedB = 20 *log10(transferFunction);
	return magnitudedB; 

}
//Metodo: Calcula la fase de 's' en grados para una frecuencia angular dada.
//--
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
    std::string _denominator;
    std::string _numerator;
    double _s_real;

    std::cout << "Introduce el numerador (por ejemplo, s^2 + 3s + 2): " << std::endl;
    std::getline(std::cin, _numerator);
    std::cout << "Introduce el denominador (por ejemplo, s^2 + 5s + 6): " << std::endl;
    std::getline(std::cin, _denominator);
    std::cout << "Introduce la frecuencia mínima (freqMin): " << std::endl;
    std::cin >> _freqMin;
    std::cout << "Introduce la frecuencia máxima (freqMax): " << std::endl;
    std::cin >> _freqMax;
    std::cout << "Introduce el valor de sigma (o) si es diferente a 0.0, de lo contrario ingresa 0.0: " << std::endl;
    std::cin >> _s_real;

    // Validación de las entradas
    if (_freqMin <= 0 || _freqMax <= _freqMin) {
        std::cerr << "Error: El rango de frecuencias no es válido. freqMin debe ser mayor a 0 y freqMax debe ser mayor a freqMin." << std::endl;
        return 1;
    }
    
    // Crear objeto de la clase MagnitudeAndPhase
    MagnitudeAndPhase mapObject(_numerator, _denominator, _freqMin, _freqMax, _s_real);

    // Procesar polos y ceros
    mapObject.processTransferFunction();
    
    // Llamar a la función frequencies()
    auto [angularFreqVector, freqVector] = mapObject.frequencies();

    // Imprimir frecuencias angulares
    std::cout << "\nFrecuencias angulares totales (w) en radianes: ";
    for (const auto& w : angularFreqVector) {
        std::cout << std::fixed << std::setprecision(4) << w << " rad "; // Agregar un espacio entre los valores
    }
    std::cout << std::endl;

    // Imprimir frecuencias en Hz
    std::cout << "Frecuencias totales (Hz): ";
    for (const auto& freq : freqVector) {
        std::cout << std::fixed << std::setprecision(4) << freq << " Hz "; // Agregar un espacio entre los valores
    }
    std::cout << std::endl;


    // Calcular e imprimir la magnitud de s para cada w
    std::cout << "\nMagnitud y Fase para cada frecuencia angular:" << std::endl;
    for (const auto& w : angularFreqVector) {
        double magnitude = mapObject.calculateMagnitude(w);
        double phase = mapObject.calculatePhase(w);

        std::cout << "  --> Para w = " << std::fixed << std::setprecision(4) << w << " rad/s:" << std::endl;
        std::cout << "      Magnitud de s = " << magnitude << std::endl;
        std::cout << "      Fase en (radianes) = " << phase << " rad" << std::endl;
        std::cout << "      Fase total en (grados) = " << phase * (180.0 / M_PI) << " grados" << std::endl;

        // Traducir el numerador y el denominador con la magnitud
        std::string translatedNumerator = mapObject.translatefunction(magnitude, true);  // Numerador
        std::string translatedDenominator = mapObject.translatefunction(magnitude, false);  // Denominador

        std::cout << "      Numerador traducido: " << translatedNumerator << std::endl;
        std::cout << "      Denominador traducido: " << translatedDenominator << std::endl;

        // Evaluar la función de transferencia
        double magnitudedB = mapObject.evaluatetranslatedfunction(magnitude);
        std::cout << "      Magnitud total en dB = " << std::fixed << std::setprecision(4) << magnitudedB << " dB" << std::endl;
    }
    std::cout << std::endl;

    return 0;
}

