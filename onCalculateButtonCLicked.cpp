#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath> // Para usar M_PI
#include "/mnt/c/Users/dravv/Projects/Qt/BodeDiagram/exprtk/exprtk.hpp"

// Función para calcular las frecuencias
std::vector<double> frecuencias(double fmin, double fmax) {
    std::vector<double> freq; 
    for (double i = fmin; i <= fmax; i++) { 
        freq.push_back(i); 
    }
    return freq;
}

// Función para calcular la magnitud en dB
double calculateMagnitude(exprtk::expression<double>& expr, double w) {
    double s_real = 0.0;
    double s_imag = w; // Representa s = j*w

// Se actualizan las variables 's_real' y 's_imag' en la tabla de símbolos
    expr.value(); // Evalúa la expresión como número complejo usando la parte imaginaria en w
    
    std::complex<double> H(s_real, s_imag);
    return 20 * std::log10(std::abs(H)); // Magnitud en dB
}
// Función para calcular la fase en grados
double calculatePhase(exprtk::expression<double>& expr, double w) {
    double s_real = 0.0;
    double s_imag = w; // Representa s = j*w

    std::complex<double> H(s_real, s_imag); 
    // Construye s como complejo usando w
    return std::arg(H) * 180.0 / M_PI; // Fase en grados
}


int main() {
    double fMin;
    double fMax;
    std::string TransferFunction;

    // Solicitar al usuario la función de transferencia
    std::cout << "Ingrese la Funcion de Trasferencia:(ej:10*s*(1000+s)/20/(20+s)/(5+s)): ";
    std::getline(std::cin, TransferFunction);

    // Solicitar al usuario los límites de frecuencia
    std::cout << "Ingrese la frecuencia mínima en (Hz): ";
    std::cin >> fMin;

    std::cout << "Ingrese la frecuencia máxima en (Hz): ";
    std::cin >> fMax;
    
    // Configuración del parser y expresión
    exprtk::symbol_table<double> symbol_table;
    symbol_table.add_constants();
    exprtk::expression<double> expr;
    expr.register_symbol_table(symbol_table);
    exprtk::parser<double> parser; 

     // Validación de la función ingresada
    if (!parser.compile(TransferFunction, expr)) {
        std::cerr << "Error en la función ingresada. Verifique la sintaxis." << std::endl;
        return 1;
    }
//Genera magnitudes
      std::vector<double> magnitudes;
      //genera fases 
    std::vector<double> phases;

    // Generar frecuencias en el rango dado
    std::vector<double> freqs = frecuencias(fMin, fMax);
    // Iterar sobre cada frecuencia
    for (double f : freqs) {
        // Calcular 'w' para cada frecuencia
        double w = 2 * M_PI * f; // w = 2 * pi * f
     double magnitude = calculateMagnitude(expr, w);
        double phase = calculatePhase(expr, w);
        // Mostrar la frecuencia y el valor completo de 'w' 
	 magnitudes.push_back(magnitude);
        phases.push_back(phase);

        std::cout << "Frecuencia f: " << f << " Hz" << std::endl;
        std::cout << "Valor de w (frecuencia angular): " << w << std::endl;
	std::cout << "Magnitud en (dB): " << magnitude <<std::endl;
	std::cout << "Fase en (grados): " << phase << std::endl; 
	
    }

    return 0;
}

