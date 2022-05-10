//
// Created by Danil-Ponkratov
//

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <electrode/electrode.hpp>

TEST_CASE( "Quick check", "[main]" ) {
  auto e = electrode(0, 2, 0, 2, -1, 3);
  int size = e.getCntAtoms();

  REQUIRE( size == 18 );
}