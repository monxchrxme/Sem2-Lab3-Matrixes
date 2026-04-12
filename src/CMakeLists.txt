# Создаем интерфейсную библиотеку matrix_lib.
# Ключевое слово INTERFACE означает, что эта библиотека header-only (состоит только из .hpp и .tpp).
# Она не компилируется в отдельный бинарник (.a или .lib), а существует только для того,
# чтобы элегантно передавать пути к заголовочным файлам всем, кто от нее зависит.
add_library(matrix_lib INTERFACE)

# Говорим CMake: "Тот, кто прилинкует matrix_lib, должен искать заголовочные файлы в текущей папке src".
target_include_directories(matrix_lib INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/lab2
)

# Создаем основной исполняемый файл нашего приложения
add_executable(Lab3_app main.cpp)

# Линкуем приложение с интерфейсной библиотекой.
# Теперь main.cpp "видит" все наши папки (memory, sequences, interfaces) без указания длинных путей.
target_link_libraries(Lab3_app PRIVATE matrix_lib)