cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(coordination-analysis VERSION 1.0.0 LANGUAGES CXX)

find_package(coretoolkit REQUIRED)

volt_add_plugin(coordination-analysis
    SOURCES src/coordination_analysis_engine.cpp
    MAIN_SOURCE src/plugin.cpp
)
