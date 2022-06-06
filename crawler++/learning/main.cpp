#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <iomanip>  // control precision
#include <memory>
#include <vector>
/*
    ofstream: write to files
    ifstream: read from files
    fstream: both read and write from/to files
*/

// int main(int argc, char const *argv[])
// {
//     std::ofstream frogFile;
//     frogFile.open("frog.txt");
//     frogFile << "hello\n" << "how are you?" << std::endl;
//     frogFile.close();
//     return 0;
// }

inline void write_something(std::ostream& os)  // can be std::cout, a file, etc.
{
    os << "hi, I am writing something to output" << std::endl;
}

int main(int argc, char const *argv[])
{
    /*
        POINTERS OVERVIEW
    */

    // REGULAR POINTERS
    int* y = nullptr;  // never set to nothing, use nullptr if needed
    y = new int[argc];
    delete[] y;

    // SMART POINTERS

    // unique pointer cannot be reassigned or converted, only _moved_ ("move semantics...")
    // unique pointer has NO OVERHEAD WHATSOEVER compared to raw pointers, and cleans itself up
    // when it goes out of scope
    // std::unique_ptr<int> x{new int}; // this initializes the unique_ptr
    // *x = 3;
    // std::cout << *x << std::endl;
    // // exit(0);

    // std::unique_ptr<double[]> froggy{new double[argc]}, froggyMoved; // initialize an array with a unique ptr
    // std::cout << froggy[argc - 1] << std::endl;
    // froggyMoved = std::move(froggy);  // transfer ownership of unique pointer to froggyMoved
    // exit(0);

    // SHARED POINTER - keeps track of references to it and cleans up once no references (a bit of overhead)
    // unlike unique_ptr, can be copied
    // std::shared_ptr<double> p1{new double};
    // std::shared_ptr<double> p3 = p1;  // this is ok
    // std::shared_ptr<double> p4 = p1;  // this is ok
    // std::cout << p4.use_count() << std::endl; // counts references
    // problems exist with cyclic references and shared_ptr, weak_ptr is a solution
    // exit(0);

    // REFERENCES
    int a = 5;
    int& b = a;  // a reference to i, kind of like a pointer, but cannot be reassigned

    // EEEEEVIL TRICK to produce undefined behavior - returning a stack allocated reference ---- oooooooh
    /*
    double& square_ref(double d)
    {
        double s = 1; // local variable, will be cleared once function returns
        return s;  // HAHA, this is eeevilll, dangling pointer
    }
    */

    // VECTORS

    const std::vector<double> vec(argc);
    std::cout << vec.size() << std::endl;
    exit(0);

    std::stringstream mysstream;
    write_something(mysstream);
    std::cout << mysstream.str() << std::endl;

    std::ifstream frogFile;
    frogFile.open("frog.txt");
    std::string contents[10];
    int i = 0;
    while (frogFile.good())
    {
        frogFile >> contents[i];   /* code */
        ++i;
    }
    // frogFile >> contents >> contents2;
    for (i = 0; i < 10; ++i)
    {
        std::cout << contents[i] << std::endl;
    }
    // std::cout << contents2 << std::endl;
    frogFile.close();

    double pi = 3.14159265359876543234567898765432;
    std::cout << std::setprecision(5) << pi << std::endl;
    return 0;
}

// // an exception type can be anything, including a struct (commonly)
// struct cannot_open_file
// {
// };

// struct lessThanZero
// {
// };


// // if exception gets thrown the program is killed if you put noexcept
// // double getSqrt(double x)
// // {
// //     double ans = -1.0;
// //     if (x > 0.0)
// //     {
// //         ans = sqrt(x);
// //     }

// //     try
// //     {
// //         if (ans == -1.0)
// //         {
// //             throw lessThanZero();
// //         }
        
// //     }
// //     catch (lessThanZero& e)
// //     {
// //         std::cout << "x must be < 0" << std::endl;
// //         exit(1);
// //     }
    
// //     return ans;
// // }

// double getSqrt(double x)
// {
//     assert(x >= 0.0);
//     return sqrt(x);
// }

// int main(int argc, char const *argv[])
// {
//     double x = atof(argv[1]);
//     x = getSqrt(x);
//     std::cout << x << std::endl;

//     bool keepTrying = true;

//     do
//     {
//         std::string fname;
//         std::cout << "please enter file name: ";
//         std::cin >> fname;
//         try 
//         {   
//             if (fname == "frog")
//             {
//                 keepTrying = false;
//             }
//             else
//             {
//                 throw cannot_open_file();
//             }
//         }
//         catch (cannot_open_file& e)
//         {
//             std::cout << "could not open file" << std::endl;
//         }
//         catch (...)
//         {
//             std::cout << "something is fishy here" << std::endl;
//         }
//     } while (keepTrying == true);

//     std::cout << "got frog!" << std::endl;
//     return 0;
// }