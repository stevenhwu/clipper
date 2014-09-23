#include <gtest/gtest.h>
//#include <gtest.h>
//#include "hippomocks.h"
#include "src/Kmer.h"
//#include "src/KmerColour.cpp"

using namespace std;


TEST(KmerColourSTest, number_of_colour){

	KmerColour colours[] =	{1,2,3,4,5,6,7,8,9};
	int expecteds[10] = 	{1,1,2,1,2,2,3,1,2};
	for (int i = 0; i < 9; ++i) {
		int no_colour = KmerColourUtil::number_of_colour_c(colours[i]);
		//	int no_colour2 = number_of_colour_s(colour);
		ASSERT_EQ(expecteds[i], no_colour);
		printf("%d\n",expecteds[i]);
	}

	KmerColour colour = 1;
	int expected = 1;
	int no_colour = KmerColourUtil::number_of_colour_c(colour);
//	int no_colour2 = number_of_colour_s(colour);
	ASSERT_EQ(expected, no_colour);
//0
	colour = 2;
	expected = 1;
	no_colour = KmerColourUtil::number_of_colour_s(colour);
	ASSERT_EQ(expected, no_colour);

	colour = 3;
	expected = 2;
	no_colour = KmerColourUtil::number_of_colour_s(colour);
	ASSERT_EQ(expected, no_colour);

	colour = 4;
	expected = 1;
	no_colour = KmerColourN::number_of_colour_n(colour);
	ASSERT_EQ(expected, no_colour);

	colour = 4;
	expected = 1;
	no_colour = KmerColourN::number_of_colour_n(colour);
	ASSERT_EQ(expected, no_colour);

	colour = 4;
	expected = 1;
	no_colour = KmerColourN::number_of_colour_n(colour);
	ASSERT_EQ(expected, no_colour);

	colour = 4;
	expected = 1;
	no_colour = KmerColourN::number_of_colour_n(colour);
	ASSERT_EQ(expected, no_colour);

	colour = 4;
	expected = 1;
	no_colour = KmerColourN::number_of_colour_n(colour);
	ASSERT_EQ(expected, no_colour);




}
//int main(){
//
//	KmerColour colour = 1;
//	int expected = 1;
//	int no_colour = KmerColourC::number_of_colour_c(colour);
//
//	colour = 4;
//	expected = 1;
//	no_colour = KmerColourC::number_of_colour_s(colour);
//	no_colour = KmerColourN::number_of_colour_n(colour);
//		//	int no_colour2 = number_of_colour_s(colour);
//	//	int no_colour2 = number_of_colour_s(colour);
//	return 0;
//}
//
//TEST(KmerColourSTest, match) {
//  ASSERT_EQ(1,1);
//  EXPECT_EQ(1,121);
//  ASSERT_EQ(1,13);
//  EXPECT_EQ(1,1231);
//}
//
//
//
//
//
//
//
//int Factorial(int n){
//	int out = 1;
//	for (int i = n; i> 0; i--) {
//		out *=i;
//	}
//	return out;// Returns the factorial of n
//}
////A test case for this function might look like:
//
//// Tests factorial of 0.
//TEST(FactorialTest, HandlesZeroInput) {
//  EXPECT_EQ(1, Factorial(0));
//}
//
//// Tests factorial of positive numbers.
//TEST(FactorialTest, HandlesPositiveInput) {
//  EXPECT_EQ(1, Factorial(1));
//  EXPECT_EQ(2, Factorial(2));
//  EXPECT_EQ(6, Factorial(3));
//  EXPECT_EQ(40320, Factorial(8));
//}


