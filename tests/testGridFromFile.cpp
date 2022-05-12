//
// Created by Danil-Ponkratov
//

#define CATCH_CONFIG_MAIN
#define bufsize 2048
#include <catch2/catch.hpp>
#include <grid/grid.hpp>

#include <cstdio>
#include <stdio.h>

TEST_CASE( "Quick check", "[main]" ) {
  auto g1 = grid(15, 10, 20, true);
  toFile("data/after_run.xyz", g1);
  auto g2 = grid(15, 10, 20, true);
  fromFile("data/after_run.xyz", g2);
  toFile("data/after_run.xyz", g2);

  char buff[bufsize];
  FILE *ch = fopen("diff data/after_run.xyz data/before_run.xyz", "r");
  if (ch == NULL)
  {
    perror("popen");
    REQUIRE(1==0);
  }
  size_t byte_count = fread(buff, 1, bufsize - 1, ch);
  printf("%s\n", buff);
  REQUIRE( byte_count == 0 );
}