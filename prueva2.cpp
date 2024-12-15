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

#define DEBUG_MODE false// need some debug? 

class MagnitudeAndPhase {
    private:
        std::string _numerator;
	std::string _denominator;
        int _freqMin;
        int _freqMax;
	double _s_real;
	std::vector<std::complex<double>> poles;
	std::vector<std::complex<double>> zeros;

	  void debug(const std::string& message) const {
        if (DEBUG_MODE) {
            std::cout << message << std::endl;
        }
    }

    public:
	bool isStable();
        MagnitudeAndPhase(std::string numerator, std::string denominator, int freqMin, int freqMax, double s_real );
	std::vector<double> parseCoefficients(const std::string& polynomial);
	std::vector<std::complex<double>> findRoots(const std::vector<double>& coefficients);
        std::pair<std::vector<double>, std::vector<double>> frequencies(); 
        std::complex<double> translateFunction(double angularFrequency, bool isNumerator);
	std::pair<double, std::complex<double>> calculateMagnitude(double angularFrequency);
        double calculatePhase(const std::complex<double>& transferFunction);
	void processTransferFunction();

};

// Constructor outside the class.
MagnitudeAndPhase::MagnitudeAndPhase(std::string numerator, std::string denominator, int freqMin, int freqMax, double s_real) 
    : _numerator(numerator), _denominator(denominator), _freqMin(freqMin), _freqMax(freqMax), _s_real(s_real){}

//--
 void MagnitudeAndPhase::processTransferFunction() {
    // Convertir cadenas a coeficientes
    debug("Debugging Method processTransferFunction: ");

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
   
//--
std::vector<double> MagnitudeAndPhase::parseCoefficients(const std::string& polynomial) {
    std::vector<double> coefficients;
    std::string modifiedPoly = polynomial;

    debug("\033[32mDebugging Method parseCoefficients: \033[0m");

    std::regex termRegex(R"(([+-]?\s*\d*\.?\d*)\s*(s(\^\d+)?)?)");
    auto termsBegin = std::sregex_iterator(modifiedPoly.begin(), modifiedPoly.end(), termRegex);
    auto termsEnd = std::sregex_iterator();

    int maxPower = 0;
    for (std::sregex_iterator i = termsBegin; i != termsEnd; ++i) {
        std::smatch match = *i;
        if (match.str().empty()) continue;

        std::string sTerm = match[2].str();
        if (!sTerm.empty()) {
            if (sTerm.find('^') != std::string::npos) {
                int power = std::stoi(sTerm.substr(sTerm.find('^') + 1));
                maxPower = std::max(maxPower, power);
            } else {
                maxPower = std::max(maxPower, 1);
            }
        }
    }

    coefficients.resize(maxPower + 1, 0.0);

    for (std::sregex_iterator i = termsBegin; i != termsEnd; ++i) {
        std::smatch match = *i;
        if (match.str().empty() || (match[1].str().empty() && match[2].str().empty())) continue;

        std::string coefStr = match[1].str();
        std::string sTerm = match[2].str();

        // Eliminar espacios en blanco
        coefStr.erase(remove_if(coefStr.begin(), coefStr.end(), ::isspace), coefStr.end());

        // Parsear coeficiente
        double coefficient = 1.0;
        if (!coefStr.empty() && coefStr != "+" && coefStr != "-") {
            coefficient = std::stod(coefStr);
        } else if (coefStr == "-") {
            coefficient = -1.0;
        }

        // Determinar potencia de s
        int power = 0;
        if (!sTerm.empty()) {
            if (sTerm.find('^') != std::string::npos) {
                power = std::stoi(sTerm.substr(sTerm.find('^') + 1));
            } else {
                power = 1;
            }
        }

        // Colocar coeficiente en la posición correcta
        coefficients[maxPower - power] = coefficient;

        // Depuración
        std::ostringstream debugMessage;
        debugMessage << "Debugging: term: " << match.str()
                     << ", coefficient: " << coefficient
                     << ", power: " << power
                     << ", maxPower: " << maxPower;
        debugMessage << "\n";
        debug(debugMessage.str());
    }

    // Depuración: Mostrar el vector final de coeficientes
    std::ostringstream debugMessage1;
    debugMessage1 << "Debugging: Coefficients vector: ";
    for (const auto& coef : coefficients) {
        debugMessage1 << coef << " ";
    }
    debugMessage1 << "\n";
    debug(debugMessage1.str());

    return coefficients;
}//--
std::vector<std::complex<double>> MagnitudeAndPhase::findRoots(const std::vector<double>& coefficients) {
    size_t degree = coefficients.size() - 1;
    debug("\033[32mDebugging Method findRoots: \033[0m");

    // Validar el grado del polinomio y los coeficientes
    if (degree < 1) throw std::invalid_argument("The polynomial degree must be at least 1.");
    if (coefficients[0] == 0) throw std::invalid_argument("The leading coefficient cannot be zero.");
    if (coefficients.size() != degree + 1) throw std::invalid_argument("The size of the coefficients vector must be degree + 1.");

    std::ostringstream debugMessage;
    debugMessage << "Debugging: Degree of polynomial: " << degree << std::endl
                 << "Debugging: Coefficients: ";
    for (const auto& coef : coefficients) {
        debugMessage << coef << " ";
    }
    debugMessage << "\n";
    debug(debugMessage.str());

    Eigen::MatrixXd companion(degree, degree);
    companion.setZero();

     
    // Asignar 1 a las subdiagonales
    for (size_t i = 1; i < degree; ++i) {
        companion( i - 1 , i ) = 1.0;
	 std::cout << "Subdiagonal: companion(" << (i - 1) << ", " << i << ") = 1.0" << std::endl;
    }
 // Asignar los coeficientes correctamente a la última fila de la matriz
    for (size_t i = 0; i < degree; ++i) {
        companion(degree - 1, i) = -coefficients[degree - i] / coefficients[0];
        std::cout << "Last row: companion(" << (degree - 1) << ", " 
                  << i << ") = " 
                  << -coefficients[degree - i] / coefficients[0] << std::endl;
    }

    debugMessage.str("");
    debugMessage << "Debugging: Companion matrix:\n" << companion << std::endl;
    debug(debugMessage.str());

    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
    Eigen::VectorXcd roots = solver.eigenvalues();

    std::vector<std::complex<double>> result(roots.size());
    for (size_t i = 0; i < roots.size(); ++i) {
        result[i] = (std::abs(_s_real) > 1e-9) ? roots[i] + _s_real : roots[i];
        debugMessage.str("");
        debugMessage << "Root before adjustment: " << roots[i] << ", after adjustment with sigma: " << result[i] << std::endl;
        debug(debugMessage.str());
    }

    debug("Debugging: Roots (final, adjusted or original):");
    for (const auto& root : result) {
        std::cout << root << std::endl;
    }

    return result;
}
//--
std::pair<std::vector<double>, std::vector<double>> MagnitudeAndPhase::frequencies() {
    std::vector<double> freqVector;
    std::vector<double> angularFreqVector;
     debug( "\033[32mDebugging Method Frequencies: \033[0m");
    //DEBUGING
    std::ostringstream debugMessage;
     debugMessage << "Debugging: Frequencies method" << std::endl
                  << "Debugging: Frequency range - Min: " << _freqMin << ", Max: " << _freqMax << std::endl;
    debug(debugMessage.str());
    //END DEBUGING
  for (double i = _freqMin; i <= _freqMax; i += 1.0) {
        freqVector.push_back(i);
	// Angular frequency calculation
        angularFreqVector.push_back(2 * M_PI * i);

        //DEBUGING
	std::ostringstream debugMessage1;
        debugMessage1 << "Frequency: " << i << ", Angular Frequency: " << angularFreqVector.back() << "rad/s."<< std::endl;
	
        debug(debugMessage1.str());
        //END DEBUGING
    }
    //DEBUGING
    std::ostringstream debugMessage2;
    debugMessage2 << "Debugging: Frequency vector size: " << freqVector.size() << std::endl
                 << "Debugging: Angular Frequency vector size: " << angularFreqVector.size() << std::endl;
    
    debug(debugMessage2.str());
    //END DEBUGING

    // Returns the pair of vectors
    return std::make_pair(angularFreqVector, freqVector); 
}
//--
std::complex<double> MagnitudeAndPhase::translateFunction(double angularFrequency, bool isNumerator) {
    std::string functionStr = isNumerator ? _numerator : _denominator;
    std::complex<double> jw(0.0, angularFrequency); // j * omega
    std::complex<double> s = (std::abs(_s_real) > 1e-9)? std::complex<double>(_s_real, angularFrequency): jw;

    // Evaluar cada término en la expresión
    std::regex termRegex("([+-]?\\d*\\.?\\d+)?(s(\\^\\d+)?)?");
    auto termsBegin = std::sregex_iterator(functionStr.begin(), functionStr.end(), termRegex);
    auto termsEnd = std::sregex_iterator();

    std::complex<double> result(0.0, 0.0); // Resultado acumulado
    debug("\033[32mDebugging Method translateFunction: \033[0m");

    std::ostringstream debugMessage;
    debugMessage << "Debugging: Evaluating function: " << functionStr << std::endl                 << "Debugging: Angular frequency (w): " << angularFrequency << std::endl;
    
    debug(debugMessage.str());

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
	std::complex<double> termValue = coefficient * std::pow(s, power);
        result += termValue;
		std::ostringstream debugMessage1;
                debugMessage1 << "Debugging: Term: " << match.str() 
                             << ", Coefficient: " << coefficient 
                             << ", Power: " << power 
                             << ", Term Value: " << termValue 
                             << ", Accumulated Result: " << result << std::endl;
                
                debug(debugMessage1.str());

    }

    std::ostringstream debugMessage2;
    debugMessage2 << "Debugging: Final translated function result: " << result << std::endl;
    
    debug(debugMessage2.str());

    return result;
}
//--
std::pair<double, std::complex<double>> MagnitudeAndPhase::calculateMagnitude(double angularFrequency) {
    // Calcular el numerador y denominador traducidos
    // Ajuste de la frecuencia compleja s = sigma + jw
    std::complex<double> s = _s_real + std::complex<double>(0, angularFrequency);

    debug("\033[32mDebugging Method CalculateMagnitude: \033[0m");

    std::ostringstream debugMessage1;
    debugMessage1 << "Debugging: s = sigma + jw = " << s << std::endl
                  << "Debugging: Magnitude of s = |s| = " << std::abs(s)
		  << std::endl;
    
    debug(debugMessage1.str());


    std::complex<double> resultNumerator = translateFunction(angularFrequency, true);
    std::complex<double> resultDenominator = translateFunction(angularFrequency, false);

    // Evitar división por un denominador cercano a cero
    if (std::abs(resultDenominator) < 1e-12) {
        std::cerr << "Error: Denominator is too close to zero!" << std::endl;
        return {};
    }
    

    // Calcular la función de transferencia (numerador / denominador)
    std::complex<double> transferFunction = resultNumerator / resultDenominator;

    // Calcular magnitud en decibeles (dB)
    double magnitudeDB = 20 * std::log10(std::abs(transferFunction));

    std::ostringstream debugMessage2;
    debugMessage2 << "Numerator: " << resultNumerator << std::endl
    		  << "Denominator: " << resultDenominator << std::endl
    		  << "Transfer Function |H(jw)| : " << transferFunction << std::endl;
                  debug(debugMessage2.str());
		   std::cout << "Transfer Function |H(jw)| : " << transferFunction 			       << std::endl
			     << "Magnitude (dB): " << magnitudeDB << std::endl;
    return {magnitudeDB, transferFunction};
}

//--
double MagnitudeAndPhase::calculatePhase(const std::complex<double>& transferFunction){
	 double phaseRadians = std::atan2(transferFunction.imag(), transferFunction.real());

    // Convertir la fase a grados
    double phaseDegrees = phaseRadians * (180.0 / M_PI);

    debug("\033[32mDebugging Method CalculatePhase: \033[0m");

    std::ostringstream debugMessage;
    debugMessage << "Debugging: Transfer Function = " << transferFunction 
	         << std::endl
                 << "Debugging: Real = " << transferFunction.real() << std::endl
                 << "Debugging: Imaginary = " << transferFunction.imag()
		 << std::endl;
                 debug(debugMessage.str());

    std::cout << "Debugging: Phase (radians) = " << phaseRadians << std::endl
              << "Debugging: Phase (degrees) = " << phaseDegrees << std::endl;

    return phaseDegrees;
}
//--
bool MagnitudeAndPhase::isStable(){
	if (poles.empty()) {
        // Si no hay polos, asumimos que el sistema no puede ser estable
        return false;
    }

for (const auto& pole : poles) { 
if (pole.real() >=0) {
return false; // unstable system if a real pole part is not negative 
}
}
return true; // stable system if all poles have a real negative part 	
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
    std::cout << "Introduce the sigma value (leave empty or enter 0 for s = jw): "              << std::endl;
    std::cin.ignore(); // Clear previous input
    std::string sigmaInput;
    std::getline(std::cin, sigmaInput);
    // Si el usuario deja vacío, usar _s_real = 0.0
    _s_real = sigmaInput.empty() ? 0.0 : std::stod(sigmaInput);

    // Validation of inputs
    if (_freqMin <= 0 || _freqMax <= _freqMin) {
        std::cerr << "Error: The frequency range is not valid. freqMin must be greater than 0 and freqMax must be greater than freqMin." << std::endl;
        return 1;
    }
    
    // Create an object of the MagnitudeAndPhase class
    MagnitudeAndPhase mapObject(_numerator, _denominator, _freqMin, _freqMax, _s_real);

    // Process poles and zeros
    mapObject.processTransferFunction();

      if (mapObject.isStable()) {
        std::cout << "The system is stable." << std::endl;
    } else {
        std::cout << "The system is unstable." << std::endl;
    }

    
  // Get frequencies and calculate magnitude and phase
    auto [angularFreqVector, freqVector] = mapObject.frequencies();
    for (const auto& w : angularFreqVector) {
        std::cout << "\nFor angular frequency w = " << std::fixed << std::setprecision(4) << w << " rad/s:\n";

	 // Calcular la magnitud
        auto[ magnitudeDB, transferFunction] = mapObject.calculateMagnitude(w);
        mapObject.calculatePhase(transferFunction);
        
    }    return 0;
}
