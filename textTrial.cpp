#include <iostream>
#include "complex.cc"
#include "input_image.cc"

int main(int argc, char* argv[])
{
	InputImage image(argv[1]);
	std::cout << image.get_width() << std::endl;
	return 0;
}
