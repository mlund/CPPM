fau_example(twobody "./" twobody.cpp)
fau_example(twobody-cuboid "./" twobody.cpp)
set_target_properties(twobody-cuboid PROPERTIES COMPILE_DEFINITIONS "CUBOID")
