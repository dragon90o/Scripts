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
    std::regex termRegex("([+-]?\\d*\\.\\d+|[+-]?\\d+)(s\\^\\d+|s)?"); 
    auto termsBegin = std::sregex_iterator(modifiedPoly.begin(), modifiedPoly.end(), termRegex);
    auto termsEnd = std::sregex_iterator();

    //  Iterate over terms and extract coefficients
    for (std::sregex_iterator i = termsBegin; i != termsEnd; ++i) {
        std::smatch match = *i;
        std::string coefStr = match[1].str();  
        std::string sTerm = match[2].str();  

        // Convert the coefficient to a number
        try {
            double coefficient = std::stod(coefStr);
            coefficients.push_back(coefficient);
        } catch (...) {
            std::cerr << "Error converting term: " << coefStr << std::endl;
        }
    }


    return coefficients;
}
    
//Method: calculates the roots (poles and zeros) of a polynomial given its coefficients.
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

    for (double i = _freqMin; i <= _freqMax; i += 1.0) {
        freqVector.push_back(i);
	// Angular frequency calculation
        angularFreqVector.push_back(2 * M_PI * i);
    }
    // Returns the pair of vectors
    return std::make_pair(angularFreqVector, freqVector); 
}
//Method: Replaces the `s` variables in the numerator or denominator with the magnitude value.
//--
std::string MagnitudeAndPhase::translatefunction(double magnitude, bool isNumerator) {
    std::string translatedfunction; 
    if (isNumerator) {
	// If processing the numerator
        translatedfunction = _numerator;  
    } else {
	// If processing the numerator
        translatedfunction = _denominator;  
    }

     // Replace all 's' with the magnitude value  
    std::string magnitudeStr = std::to_string(magnitude);
    size_t pos = 0;
    while ((pos = translatedfunction.find("s", pos)) != std::string::npos) {
        translatedfunction.replace(pos, 1, "(" + magnitudeStr + ")");
        pos += magnitudeStr.length();
    }

    return translatedfunction;
}

//Method: Calculates the magnitude s, "|s| = |j * w|" for a given angular frequency.
//--
double MagnitudeAndPhase::calculateMagnitude(double w) {
    
    std::complex<double> s(_s_real, w);  // s = j * w
    return std::abs(s);  // Retorna la magnitud |s|
}

//Method: Evaluates the transfer function with the given magnitude (numerator/denominator).
//--
double MagnitudeAndPhase::evaluatetranslatedfunction(double magnitude){

        std::string translatedNumerator = translatefunction(magnitude, true);
        std::string translatedDenominator = translatefunction(magnitude, false);

	// Evaluate numerator
	exprtk::expression<double> expressionNumerator;
	exprtk::parser<double> parserNumerator;
	if (!parserNumerator.compile(translatedNumerator, expressionNumerator)) {
		std::cerr << "Error compiling the numerator: "<<parserNumerator.error() << std::endl;
		return 0.0; 
	}
	double resultNumerator = expressionNumerator.value();

	//Evaluate denominator
	exprtk::expression<double> expressionDenominator;
	exprtk::parser<double> parserDenominator;
	if (!parserDenominator.compile(translatedDenominator, expressionDenominator)) {
        std::cerr << "Error compiling the denominator: " << parserDenominator.error() << std::endl;
        return 0.0;
    }
	double resultDenominator = expressionDenominator.value();

	double transferFunction= resultNumerator / resultDenominator;

	double magnitudedB = 20 *log10(transferFunction);
	return magnitudedB; 

}
//Method: Calculates the phase of 's' in degrees for a given angular frequency.
//--
double MagnitudeAndPhase::calculatePhase(double w){
	std::complex<double> s(_s_real, w);
        double magnitude = std::abs(s);
        double phase = std::arg(s);

        std::cout << "s: " << s << " -> Magnitude: " << magnitude << ", Phase (radians)?: " << phase << std::endl;

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
    
    // Call the frequencies() function
    auto [angularFreqVector, freqVector] = mapObject.frequencies();

    // Print angular frequencies
    std::cout << "\nTotal angular frequencies (w) in radians: ";
    for (const auto& w : angularFreqVector) {
        std::cout << std::fixed << std::setprecision(4) << w << " rad "; // Agregar un espacio entre los valores
    }
    std::cout << std::endl;

    // Print frequencies in Hz
    std::cout << "Total frequencies (Hz): ";
    for (const auto& freq : freqVector) {
        std::cout << std::fixed << std::setprecision(4) << freq << " Hz "; // Agregar un espacio entre los valores
    }
    std::cout << std::endl;


    // Calculate and print the magnitude of s for each w
    std::cout << "\nMagnitude and Phase for each angular frequency: " << std::endl;
    for (const auto& w : angularFreqVector) {
        double magnitude = mapObject.calculateMagnitude(w);
        double phase = mapObject.calculatePhase(w);

        std::cout << "  --> For w = " << std::fixed << std::setprecision(4) << w << " rad/s:" << std::endl;
        std::cout << "      Magnitude of s = " << magnitude << std::endl;
        std::cout << "      Phase (in degrees) = " << phase << " degrees" << std::endl;
        

        // Translate the numerator and denominator with the magnitude
	// Numerator
	std::string translatedNumerator = mapObject.translatefunction(magnitude, true);
	//Denominator
        std::string translatedDenominator = mapObject.translatefunction(magnitude, false); 

        std::cout << "      Translated numerator: " << translatedNumerator << std::endl;
        std::cout << "       Translated denominator: " << translatedDenominator << std::endl;

        // Evaluate the transfer function
	double magnitudedB = mapObject.evaluatetranslatedfunction(magnitude);
        std::cout << "      Total magnitude in dB = " << std::fixed << std::setprecision(4) << magnitudedB << " dB" << std::endl;
    }
    std::cout << std::endl;

    return 0;
}

