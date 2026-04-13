# Создаем исполняемый файл для тестов. Сюда мы передаем наш файл с макросами TEST()
add_executable(run_tests
        test_complex.cpp
        test_matrix.cpp
        test_square_matrix.cpp
        test_vector.cpp
        test_diagonal_matrix.cpp
        test_triangular_matrix.cpp
        test_system_of_equations.cpp
)

# Самая важная магия линковки:
# 1. Прикрепляем нашу библиотеку matrix_lib, чтобы тесты "видели" классы DynamicArray, ListSequence и т.д.
# 2. Прикрепляем GTest::gtest — саму библиотеку Google Test (макросы EXPECT_EQ, EXPECT_THROW).
# 3. Прикрепляем GTest::gtest_main — это готовая реализация функции int main(),
#    которая сама найдет все тесты и запустит их.
target_link_libraries(run_tests PRIVATE
        matrix_lib
        GTest::gtest
        GTest::gtest_main
)
# Сообщаем CTest, что бинарник run_tests содержит тесты Google Test.
# Это позволит IDE (CLion, VS Code, Visual Studio) красиво отображать список тестов
# в графическом интерфейсе, с зелеными галочками при успехе.
gtest_discover_tests(run_tests)