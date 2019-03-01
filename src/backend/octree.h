#pragma once

#include <limits>
#include <cmath>
#include <vector>
#include <tuple>

namespace inspire {
  // NOTE: It could also be in common interface, but as it is tailor-made for molecules, it is placed here
  namespace backend {
    class OctreeNode {
      protected:
      static const size_t LIMIT = 4096;

      inline double squared(const double d) { return d*d; }

      public:
      virtual bool contact(const std::tuple<double, double, double> &point, const double distance) = 0;
    };

    class OctreeList : public OctreeNode {
      private:
      const std::vector<std::pair<std::tuple<double, double, double>, double>> POINTS;

      public:
      OctreeList(const std::vector<std::pair<std::tuple<double, double, double>, double>> &points) : POINTS(points) { }

      ~OctreeList() {
        // TODO: Check whether is more effective this or const POINTS
        // POINTS.clear();
      }

      bool contact(const std::tuple<double, double, double> &point, const double distance) override {
        for (auto points_it = POINTS.begin(); points_it != POINTS.end(); ++points_it) {
          if (squared(std::get<0>(point)-std::get<0>(points_it->first)) + squared(std::get<1>(point)-std::get<1>(points_it->first)) +
              squared(std::get<2>(point)-std::get<2>(points_it->first)) <= squared(distance+points_it->second)) {
            return true;
          }
        }
        return false;
      }
    };

    class OctreeInnerNode : public OctreeNode {
      private:
      const std::tuple<double, double, double> SPLITTING_POINT;
      std::pair<OctreeNode*, double> NODES[2][2][2];

      public:
      OctreeInnerNode(const std::vector<std::pair<std::tuple<double, double, double>, double>> &points, const double x, const double y, const double z) : SPLITTING_POINT(x, y, z) {
        std::vector<std::pair<std::tuple<double, double, double>, double>> nodes[2][2][2];
        std::tuple<std::pair<double, double>, std::pair<double, double>, std::pair<double, double>> boundaries[2][2][2];
        for (size_t i = 0; i < 2; i++) {
          for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
              NODES[i][j][k].second = 0;
              std::get<0>(boundaries[i][j][k]).first = std::numeric_limits<double>::max();
              std::get<0>(boundaries[i][j][k]).second = std::numeric_limits<double>::lowest();
              std::get<1>(boundaries[i][j][k]).first = std::numeric_limits<double>::max();
              std::get<1>(boundaries[i][j][k]).second = std::numeric_limits<double>::lowest();
              std::get<2>(boundaries[i][j][k]).first = std::numeric_limits<double>::max();
              std::get<2>(boundaries[i][j][k]).second = std::numeric_limits<double>::lowest();
            }
          }
        }

        for (auto points_it = points.begin(); points_it != points.end(); ++points_it) {
          size_t i = std::get<0>(points_it->first) < std::get<0>(SPLITTING_POINT) ? 0 : 1;
          size_t j = std::get<1>(points_it->first) < std::get<1>(SPLITTING_POINT) ? 0 : 1;
          size_t k = std::get<2>(points_it->first) < std::get<2>(SPLITTING_POINT) ? 0 : 1;
          nodes[i][j][k].push_back(*points_it);
          if (std::get<0>(boundaries[i][j][k]).first > std::get<0>(points_it->first)) {
            std::get<0>(boundaries[i][j][k]).first = std::get<0>(points_it->first);
          }
          if (std::get<0>(boundaries[i][j][k]).second < std::get<0>(points_it->first)) {
            std::get<0>(boundaries[i][j][k]).second = std::get<0>(points_it->first);
          }
          if (std::get<1>(boundaries[i][j][k]).first > std::get<1>(points_it->first)) {
            std::get<1>(boundaries[i][j][k]).first = std::get<1>(points_it->first);
          }
          if (std::get<1>(boundaries[i][j][k]).second < std::get<1>(points_it->first)) {
            std::get<1>(boundaries[i][j][k]).second = std::get<1>(points_it->first);
          }
          if (std::get<2>(boundaries[i][j][k]).first > std::get<2>(points_it->first)) {
            std::get<2>(boundaries[i][j][k]).first = std::get<2>(points_it->first);
          }
          if (std::get<2>(boundaries[i][j][k]).second < std::get<2>(points_it->first)) {
            std::get<2>(boundaries[i][j][k]).second = std::get<2>(points_it->first);
          }
          if (points_it->second > NODES[i][j][k].second) {
            NODES[i][j][k].second = points_it->second;
          }
        }

        for (size_t i = 0; i < 2; i++) {
          for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
              if (nodes[i][j][k].size() > OctreeNode::LIMIT && (std::get<0>(boundaries[i][j][k]).first != std::get<0>(boundaries[i][j][k]).second ||
                                                                std::get<1>(boundaries[i][j][k]).first != std::get<1>(boundaries[i][j][k]).second ||
                                                                std::get<2>(boundaries[i][j][k]).first != std::get<2>(boundaries[i][j][k]).second)) {
                NODES[i][j][k].first = new OctreeInnerNode(nodes[i][j][k], (std::get<0>(boundaries[i][j][k]).first + std::get<0>(boundaries[i][j][k]).second)/2,
                                                                           (std::get<1>(boundaries[i][j][k]).first + std::get<1>(boundaries[i][j][k]).second)/2,
                                                                           (std::get<2>(boundaries[i][j][k]).first + std::get<2>(boundaries[i][j][k]).second)/2);
              } else {
                NODES[i][j][k].first = new OctreeList(nodes[i][j][k]);
              }
            }
          }
        }
      }

      ~OctreeInnerNode() {
        for (size_t i = 0; i < 2; i++) {
          for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
              if (NODES[i][j][k].first != nullptr) {
                delete NODES[i][j][k].first;
                NODES[i][j][k].first = nullptr;
              }
            }
          }
        }
      }

      bool contact(const std::tuple<double, double, double> &point, const double distance) override {
        const size_t x = std::get<0>(point) < std::get<0>(SPLITTING_POINT) ? 0 : 1;
        const size_t y = std::get<1>(point) < std::get<1>(SPLITTING_POINT) ? 0 : 1;
        const size_t z = std::get<2>(point) < std::get<2>(SPLITTING_POINT) ? 0 : 1;
        // This order was choosen because of probability of collisions;
        // TODO: Check efficiency of varying order of axis based on the size of intersection
        // The same cube
        if (NODES[x][y][z].first->contact(point, distance)) {
          return true;
        }
        // Side-neighboring cubes
        if (std::abs(std::get<0>(point) - std::get<0>(SPLITTING_POINT)) <= distance + NODES[1-x][y][z].second && NODES[1-x][y][z].first->contact(point, distance)) {
          return true;
        }
        if (std::abs(std::get<1>(point) - std::get<1>(SPLITTING_POINT)) <= distance + NODES[x][1-y][z].second && NODES[x][1-y][z].first->contact(point, distance)) {
          return true;
        }
        if (std::abs(std::get<2>(point) - std::get<2>(SPLITTING_POINT)) <= distance + NODES[x][y][1-z].second && NODES[x][y][1-z].first->contact(point, distance)) {
          return true;
        }
        // Edge-neighboring cubes
        if (squared(std::get<0>(point)-std::get<0>(SPLITTING_POINT)) + squared(std::get<1>(point)-std::get<1>(SPLITTING_POINT)) <= squared(distance + NODES[1-x][1-y][z].second) &&
            NODES[1-x][1-y][z].first->contact(point, distance)) {
          return true;
        }
        if (squared(std::get<0>(point)-std::get<0>(SPLITTING_POINT)) + squared(std::get<2>(point)-std::get<2>(SPLITTING_POINT)) <= squared(distance + NODES[1-x][y][1-z].second) &&
            NODES[1-x][y][1-z].first->contact(point, distance)) {
          return true;
        }
        if (squared(std::get<1>(point)-std::get<1>(SPLITTING_POINT)) + squared(std::get<2>(point)-std::get<2>(SPLITTING_POINT)) <= squared(distance + NODES[x][1-y][1-z].second) &&
            NODES[x][1-y][1-z].first->contact(point, distance)) {
          return true;
        }
        // Vertex-neighboring cube
        if (squared(std::get<0>(point)-std::get<0>(SPLITTING_POINT)) + squared(std::get<1>(point)-std::get<1>(SPLITTING_POINT)) + squared(std::get<2>(point)-std::get<2>(SPLITTING_POINT)) <=
            squared(distance + NODES[1-x][1-y][1-z].second) && NODES[1-x][1-y][1-z].first->contact(point, distance)) {
          return true;
        }

        return false;
      }
    };

    class Octree : public OctreeNode {
      private:
      // Definition of cuboid with points [{min_x, max_x}, {min_y, max_y}, {min_z, max_z}]
      std::tuple<std::pair<double, double>, std::pair<double, double>, std::pair<double, double>> BOUNDARIES;
      // Max radius of points in octree
      double DISTANCE;
      // Root of octree
      OctreeNode* ROOT;

      public:
      Octree(const std::vector<std::pair<std::tuple<double, double, double>, double>> &points) : BOUNDARIES({std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()},
                                                                                                     {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()},
                                                                                                     {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()}), DISTANCE(0) {
        for (auto points_it = points.begin(); points_it != points.end(); ++points_it) {
          if (std::get<0>(points_it->first) < std::get<0>(BOUNDARIES).first) {
            std::get<0>(BOUNDARIES).first = std::get<0>(points_it->first);
          }
          if (std::get<0>(points_it->first) > std::get<0>(BOUNDARIES).second) {
            std::get<0>(BOUNDARIES).second = std::get<0>(points_it->first);
          }
          if (std::get<1>(points_it->first) < std::get<1>(BOUNDARIES).first) {
            std::get<1>(BOUNDARIES).first = std::get<1>(points_it->first);
          }
          if (std::get<1>(points_it->first) > std::get<1>(BOUNDARIES).second) {
            std::get<1>(BOUNDARIES).second = std::get<1>(points_it->first);
          }
          if (std::get<2>(points_it->first) < std::get<2>(BOUNDARIES).first) {
            std::get<2>(BOUNDARIES).first = std::get<2>(points_it->first);
          }
          if (std::get<2>(points_it->first) > std::get<2>(BOUNDARIES).second) {
            std::get<2>(BOUNDARIES).second = std::get<2>(points_it->first);
          }
          if (points_it->second > DISTANCE) {
            DISTANCE = points_it->second;
          }
        }

        if (points.size() > LIMIT && (std::get<0>(BOUNDARIES).first != std::get<0>(BOUNDARIES).second || std::get<1>(BOUNDARIES).first != std::get<1>(BOUNDARIES).second ||
                                      std::get<2>(BOUNDARIES).first != std::get<2>(BOUNDARIES).second)) {
          ROOT = new OctreeInnerNode(points, (std::get<0>(BOUNDARIES).first + std::get<0>(BOUNDARIES).second) / 2,
                                             (std::get<1>(BOUNDARIES).first + std::get<1>(BOUNDARIES).second) / 2,
                                             (std::get<2>(BOUNDARIES).first + std::get<2>(BOUNDARIES).second) / 2);
        } else {
          ROOT = new OctreeList(points);
        }
      }

      ~Octree() {
        if (ROOT != nullptr) {
          delete ROOT;
          ROOT = nullptr;
        }
      }

      bool contact(const std::tuple<double, double, double> &point, const double distance) override {
        if (std::get<0>(point)+distance < std::get<0>(BOUNDARIES).first-DISTANCE) {
          return false;
        }
        if (std::get<0>(BOUNDARIES).second+DISTANCE < std::get<0>(point)-distance) {
          return false;
        }
        if (std::get<1>(point) + distance < std::get<1>(BOUNDARIES).first - DISTANCE) {
          return false;
        }
        if (std::get<1>(BOUNDARIES).second + DISTANCE < std::get<1>(point) - distance) {
          return false;
        }
        if (std::get<2>(point) + distance < std::get<2>(BOUNDARIES).first - DISTANCE) {
          return false;
        }
        if (std::get<2>(BOUNDARIES).second + DISTANCE < std::get<2>(point) - distance) {
          return false;
        }
        // TODO: It should not be efficient to check diagonals too, but test it
        return ROOT->contact(point, distance);
      }
    };
  }
}
