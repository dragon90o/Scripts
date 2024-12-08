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
        std::complex<double> translateFunction(double angularFrequency, bool isNumerator);
	void calculateMagnitude(double angularFrequency);
        double calculatePhase(double w);
	void processTransferFunction();

};

// Constructor outside the class.
MagnitudeAndPhase::MagnitudeAndPhase(std::string numerator, std::string denominator, int freqMin, int freqMax, double s_real) 
    : _numerator(numerator), _denominator(denominator), _freqMin(freqMin), _freqMax(freqMax), _s_real(s_real){}

//Method: converts a polynomial string into a list of numerical coefficients. 
//-- polynomial is the instance to pass numerator and denominator
std::vector<double> MagnitudeAndPhase::parseCoefficients(const std::string& polynomial) {
    std::vector<double> coefficients;
    std::string modifiedPoly = polynomial;

    
     // Use a regular expression to separate the expression and properly handle the operators.
     // Regex to handle coefficients with s^n or s
    std::regex termRegex("([+-]?\\d*\\.?\\d+)?(s(\\^\\d+)?)?");
    auto termsBegin = std::sregex_iterator(modifiedPoly.begin(), modifiedPoly.end(), termRegex);
    auto termsEnd = std::sregex_iterator();
//DEBUGING -->
     std::cout << "Debugging: Found " << std::distance(termsBegin, termsEnd) << " matches." << std::endl;
//END DEBBUGING -->
    //  Iterate over terms and extract coefficients
    int maxPower = 0;
    for (std::sregex_iterator i = termsBegin; i != termsEnd; ++i) {
        std::smatch match = *i;
	if (match.str().empty()) continue;

        std::string sTerm = match[2].str();  

        if (!sTerm.empty()) { // Check if it's a term with 's'
            if (sTerm.find('^') != std::string::npos) {
                int power = std::stoi(sTerm.substr(sTerm.find('^') + 1));
                maxPower = std::max(maxPower, power);
            } else {
                maxPower = std::max(maxPower, 1); // 's' without '^n' means power 1
            }
        }
    }

    // Initialize the coefficients vector with zeros
    coefficients.resize(maxPower + 1, 0.0);

    // Fill coefficients in the proper positions
    for (std::sregex_iterator i = termsBegin; i != termsEnd; ++i) {
        std::smatch match = *i;
	if (match.str().empty() || (match[1].str().empty() && match[2].str().empty())) continue;

        std::string coefStr = match[1].str();  // Coefficient part
        std::string sTerm = match[2].str();    // 's' term

        // Parse coefficient
        double coefficient = 1.0; // Default is 1
        if (!coefStr.empty()) {
            coefficient = std::stod(coefStr);
        }

        // Determine power of s
        int power = 0;
        if (!sTerm.empty()) {
            if (sTerm.find('^') != std::string::npos) {
                power = std::stoi(sTerm.substr(sTerm.find('^') + 1));
            } else {
                power = 1; // 's' without '^n' means power 1
            }
        }

        // Place coefficient in the correct position
        coefficients[maxPower - power] = coefficient;
//DEBUGING
        // Debugging: Show each term
        std::cout << "Debugging: term: " << match.str()
                  << ", coefficient: " << coefficient
                  << ", power: " << power
                  << ", maxPower: " << maxPower << std::endl;
    }

    // Debugging: Show the final coefficients vector
    std::cout << "Debugging: Coefficients vector: ";
    for (const auto& coef : coefficients) {
        std::cout << coef << " ";
    }
    std::cout << std::endl;
//END DEBUGING
    return coefficients;
}    
//Method: calculates the roots (poles and zeros) of a polynomial given its coefficients.
//--
std::vector<std::complex<double>> MagnitudeAndPhase::findRoots(const std::vector<double>& coefficients) {
    int degree = coefficients.size() - 1;
    //DEBUGING -->
    std::cout << "DEBUGING METODO FINDROOTS"<< std::endl;
    std::cout << "Debugging: Degree of polynomial: " << degree << std::endl;
    // Mostrar los coeficientes
    std::cout << "Debugging: Coefficients: ";
    for (const auto& coef : coefficients) {
        std::cout << coef << " ";
    }
    std::cout << std::endl;

    // END DEBUGING -->

    Eigen::MatrixXd companion(degree, degree);
    companion.setZero();

for (int i = 0; i < degree; ++i) {
        // Asignar coeficientes en la última columna
        std::cout << "Assigning coefficient[" << degree - i << "] = "
                  << coefficients[degree - i] << " to companion matrix.\n";

        companion(i, degree - 1) = -coefficients[degree - i] / coefficients.front();

        // Asignar las subdiagonales
        if (i < degree - 1) {
            companion(i + 1, i) = 1.0;
        }
    }
//DEBUGING --> 
   // Mostrar la matriz compañera
    std::cout << "Debugging: Companion matrix:\n" << companion << std::endl;
//END DEBUGING-->
    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
    Eigen::VectorXcd roots = solver.eigenvalues();
//DEBUGING -->
 // Mostrar raíces calculadas
    std::cout << "Debugging: Roots (eigenvalues):" << std::endl;
    for (int i = 0; i < roots.size(); ++i) {
        std::cout << roots[i] << std::endl;
    }
//END DEBUGING -->
    std::vector<std::complex<double>> result(roots.size());
    for (int i = 0; i < roots.size(); ++i) {
        result[i] = roots[i];
    }
    return result;
}
//Method: Processes the transfer function: calculates zeros, poles, and prints the results.
//--
void MagnitudeAndPhase::processTransferFunction() {
    // Convertir cadenas a coeficientes
    auto numCoefficients = parseCoefficients(_numerator);
    auto denCoefficients = parseCoefficients(_denominator);

    // Calculate zeros and poles
    zeros = findRoots(numCoefficients);
    poles = findRoots(denCoefficients);

    //Print results
    std::cout << "Zeros (Roots of the Numerator): ";
    for (const auto& zero : zeros) {
        std::cout << zero << " ";
    }
    std::cout << std::endl;

    std::cout << "Poles (Roots of the Denominator): ";
    for (const auto& pole : poles) {
        std::cout << pole << " ";
    }
    std::cout << std::endl;
}

//Method: generates two vectors; frequencies in Hz and angular frequencies (w).
//--
std::pair<std::vector<double>, std::vector<double>> MagnitudeAndPhase::frequencies() {
    std::vector<double> freqVector;
    std::vector<double> angularFreqVector;
//DEBUGING
    std::cout << "DEBUGING FREQUENCIES"<< std::endl;
     // Debug: Verificar límites
    std::cout << "Debugging: Frequencies method" << std::endl;
    std::cout << "Debugging: Frequency range - Min: " << _freqMin << ", Max: " << _freqMax << std::endl;
//END DEBUGING
    for (double i = _freqMin; i <= _freqMax; i += 1.0) {
        freqVector.push_back(i);
	// Angular frequency calculation
        angularFreqVector.push_back(2 * M_PI * i);
//DEBUGING
	 // Debug: Mostrar valores en cada iteración
        std::cout << "Frequency: " << i << ", Angular Frequency: " << angularFreqVector.back() << "rad/s."<< std::endl;
//END DEBUGING
    }
    //DEBUGING
    // Debug: Mostrar tamaños de los vectores
    std::cout << "Debugging: Frequency vector size: " << freqVector.size() << std::endl;
    std::cout << "Debugging: Angular Frequency vector size: " << angularFreqVector.size() << std::endl;
//END DEBUGING
    // Returns the pair of vectors
    return std::make_pair(angularFreqVector, freqVector); 
}
//Method: Replaces the `s` variables in the numerator or denominator with the magnitude value.
//--
std::complex<double> MagnitudeAndPhase::translateFunction(double angularFrequency, bool isNumerator) {
    std::string functionStr = isNumerator ? _numerator : _denominator;
    std::complex<double> jw(0.0, angularFrequency); // j * omega

    // Evaluar cada término en la expresión
    std::regex termRegex("([+-]?\\d*\\.?\\d+)?(s(\\^\\d+)?)?");
    auto termsBegin = std::sregex_iterator(functionStr.begin(), functionStr.end(), termRegex);
    auto termsEnd = std::sregex_iterator();

    std::complex<double> result(0.0, 0.0); // Resultado acumulado

    // Mostrar la función a evaluar
    std::cout << "Debugging: Evaluating function: " << functionStr << std::endl;
    std::cout << "Debugging: Angular frequency (w): " << angularFrequency << std::endl;

    for (auto it = termsBegin; it != termsEnd; ++it) {
        std::smatch match = *it;
        if (match.str().empty()) continue;

        // Obtener coeficiente
        double coefficient = 1.0; // Por defecto
        if (!match[1].str().empty()) {
            coefficient = std::stod(match[1].str());
        }

        // Obtener potencia de s
        int power = 0;
        if (!match[2].str().empty()) {
            if (match[2].str().find('^') != std::string::npos) {
                power = std::stoi(match[2].str().substr(match[2].str().find('^') + 1));
            } else {
                power = 1;
            }
        }

        // Evaluar el término
        std::complex<double> termValue = coefficient * std::pow(jw, power);
        result += termValue;

        // Mostrar información de depuración por cada término
        std::cout << "Debugging: Term: " << match.str() 
                  << ", Coefficient: " << coefficient 
                  << ", Power: " << power 
                  << ", Term Value: " << termValue 
                  << ", Accumulated Result: " << result 
                  << std::endl;
    }

    // Mostrar resultado final
    std::cout << "Debugging: Final translated function result: " << result << std::endl;

    return result;
}
  void MagnitudeAndPhase::calculateMagnitude(double angularFrequency) {
    // Calcular el numerador y denominador traducidos
    std::complex<double> resultNumerator = translateFunction(angularFrequency, true);
    std::complex<double> resultDenominator = translateFunction(angularFrequency, false);

    // Evitar división por un denominador cercano a cero
    if (std::abs(resultDenominator) < 1e-12) {
        std::cerr << "Error: Denominator is too close to zero!" << std::endl;
        return;
    }

    // Calcular la función de transferencia (numerador / denominador)
    std::complex<double> transferFunction = resultNumerator / resultDenominator;

    // Calcular magnitud en decibeles (dB)
    double magnitudeDB = 20 * std::log10(std::abs(transferFunction));

    // Mostrar resultados
    std::cout << "Numerator: " << resultNumerator << std::endl;
    std::cout << "Denominator: " << resultDenominator << std::endl;
    std::cout << "Transfer Function: " << transferFunction << std::endl;
    std::cout << "Magnitude (dB): " << magnitudeDB << std::endl;
}



//Method: Calculates the phase of 's' in degrees for a given angular frequency.
//--
double MagnitudeAndPhase::calculatePhase(double w){
	std::complex<double> s(_s_real, w);
        double angularFrequency = std::abs(s);
        double phase = std::arg(s);

        std::cout << "s: " << s << " -> angularFrequency: " << angularFrequency  << ", Phase (radians)?: " << phase << std::endl;

	// Returns the phase in degrees
	return std::arg(s)*(180.0 /M_PI);
}


int main() {
    double _freqMin;
    double _freqMax;
    std::string _denominator;
    std::string _numerator;
    double _s_real;

    std::cout << "Introduce the numerator (e.g., s^2 + 3s + 2): " << std::endl;
    std::getline(std::cin, _numerator);
    std::cout << "Introduce the denominator (e.g., s^2 + 5s + 6): " << std::endl;
    std::getline(std::cin, _denominator);
    std::cout << "Introduce the minimum frequency (freqMin): " << std::endl;
    std::cin >> _freqMin;
    std::cout << "Introduce the maximum frequency (freqMax): " << std::endl;
    std::cin >> _freqMax;
    std::cout << "Introduce the sigma value (or) if it is different from 0.0, otherwise enter 0: " << std::endl;
    std::cin >> _s_real;

    // Validation of inputs
    if (_freqMin <= 0 || _freqMax <= _freqMin) {
        std::cerr << "Error: The frequency range is not valid. freqMin must be greater than 0 and freqMax must be greater than freqMin." << std::endl;
        return 1;
    }
    
    // Create an object of the MagnitudeAndPhase class
    MagnitudeAndPhase mapObject(_numerator, _denominator, _freqMin, _freqMax, _s_real);

    // Process poles and zeros
    mapObject.processTransferFunction();
    
  // Get frequencies and calculate magnitude and phase
    auto [angularFreqVector, freqVector] = mapObject.frequencies();
    for (const auto& w : angularFreqVector) {
        std::cout << "\nFor angular frequency w = " << std::fixed << std::setprecision(4) << w << " rad/s:\n";
        mapObject.calculateMagnitude(w);
        double phase = mapObject.calculatePhase(w);
        std::cout << "  Phase (degrees): " << phase << "\n";
    }    return 0;
}

