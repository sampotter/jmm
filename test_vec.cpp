#include <catch2/catch.hpp>

#include "vec.h"

using namespace Catch::literals;

TEST_CASE ("dvec3_dbl_div works", "[dvec3]") {
  dvec3 u = {.data = {1.0, 2.0, 3.0}};
  dbl a = 0.5;
  dvec3 v = dvec3_dbl_div(u, a);
  REQUIRE(v.data[0] == Approx(2.0));
  REQUIRE(v.data[1] == Approx(4.0));
  REQUIRE(v.data[2] == Approx(6.0));
}

TEST_CASE ("dvec3_dist works", "[dvec3]") {
  dvec3 u = {
    .data = {
      0.6417927091910892,
      -1.1362528112059935,
      -0.8477659762445305
    }
  };
  dvec3 v = {
    .data = {
      -1.0566292112977629,
      0.81411758780935,
      0.17566143340868834
    }
  };
  dbl d = dvec3_dist(u, v);
  REQUIRE(d == Approx(2.781363941698714));
}

TEST_CASE ("dvec3_dot works", "[dvec3]") {
  dvec3 u = {
    .data = {
      -0.16843202882104324,
      -2.1290157800245857,
      -3.0639551665912337
    }
  };
  dvec3 v = {
    .data = {
      -0.5885364261856167,
      -0.5736500179311756,
      0.16633253999640382
    }
  };
  REQUIRE(dvec3_dot(u, v) == Approx(0.8108028793901646));
}

TEST_CASE ("dvec3_maxdist works", "[dvec3]") {
  dvec3 u = {
    .data = {
      -1.1049973403391797,
      1.0480109216157314,
      -0.3075628963377684
    }
  };
  dvec3 v = {
    .data = {
      -0.08421444863196027,
      1.0322275154440248,
      0.9109219010584781
    }
  };
  REQUIRE(dvec3_maxdist(u, v) == Approx(1.2184847973962465));
}

TEST_CASE ("dvec3_maxnorm works", "[dvec3]") {
  dvec3 u = {
    .data = {
      0.0013148458299276043,
      -1.2592122006954491,
      1.0648539656693685
    }
  };
  REQUIRE(dvec3_maxnorm(u) == Approx(1.2592122006954491));
}

TEST_CASE ("dvec3_nan works", "[dvec3]") {
  dvec3 nan = dvec3_nan();
  REQUIRE(isnan(nan.data[0]));
  REQUIRE(isnan(nan.data[1]));
  REQUIRE(isnan(nan.data[2]));
}

TEST_CASE ("dvec3_norm works", "[dvec3]") {
  dvec3 u = {
    .data = {
      -0.434521508703326,
      -0.019900125495633626,
      -1.6719997244423974
    }
  };
  REQUIRE(dvec3_norm(u) == Approx(1.7276539106707713));
}

TEST_CASE ("dvec3_normalized works", "[dvec3]") {
  dvec3 u = {
    .data = {
      0.4206872480720245,
      0.6121063429157634,
      -0.5368203638129566
    }
  };
  u = dvec3_normalized(u);
  REQUIRE(u.data[0] == Approx(0.45905440886642385));
  REQUIRE(u.data[1] == Approx(0.6679311452827229));
  REQUIRE(u.data[2] == Approx(-0.5857789983104621));
}

TEST_CASE ("dvec3_saxpy works", "[dvec3]") {
  dvec3 u = {
    .data = {
      -0.35850101736644807,
      -0.22715015123286572,
      0.7616333022658741
    }
  };
  dvec3 v = {
    .data = {
      -0.1776469197740194,
      -1.2767905968125333,
      -1.8875953611972274
    }
  };
  dbl a = 1.4274190019383162;
  dvec3 w = dvec3_saxpy(a, u, v);
  REQUIRE(w.data[0] == Approx(-0.6893780841771057));
  REQUIRE(w.data[1] == Approx(-1.601029038975488));
  REQUIRE(w.data[2] == Approx(-0.8004255130338895));
}

TEST_CASE ("dvec3_sub works", "[dvec3]") {
  dvec3 u = {
    .data = {
      -0.8266017380400008,
      1.8406146995450003,
      0.3377554584592156
    }
  };
  dvec3 v = {
    .data = {
      -1.6890253296151674,
      1.243353922818447,
      0.5372978011342767
    }
  };
  dvec3 w = dvec3_sub(u, v);
  REQUIRE(w.data[0] == Approx(0.8624235915751666));
  REQUIRE(w.data[1] == Approx(0.5972607767265534));
  REQUIRE(w.data[2] == Approx(-0.1995423426750611));
}

TEST_CASE ("ivec3_dbl_mul works", "[ivec3]") {
  ivec3 p = {.data = {1, 2, 3}};
  dbl a = {0.5};
  dvec3 u = ivec3_dbl_mul(p, a);
  REQUIRE(u.data[0] == Approx(0.5));
  REQUIRE(u.data[1] == Approx(1.0));
  REQUIRE(u.data[2] == Approx(1.5));
}
