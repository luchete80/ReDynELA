#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

// A helper function to create a sorted pair for an edge
std::pair<int, int> createEdge(int a, int b) {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}


int test1() {
    // Define the mesh as a vector of triangles, each represented by 3 vertex indices
    std::vector<std::vector<int>> triangles = {
        {1, 2, 3},
        {3, 2, 4},
        {4, 2, 5},
        {5, 2, 1}
    };

    // Map to store edges and their counts
    std::map<std::pair<int, int>, int> edgeCount;

    // Iterate over each triangle and its edges
    for (const auto& triangle : triangles) {
        for (int i = 0; i < 3; ++i) {
            int v1 = triangle[i];
            int v2 = triangle[(i + 1) % 3]; // Next vertex in the triangle
            std::pair<int, int> edge = createEdge(v1, v2);
            edgeCount[edge]++;
        }
    }

    // Print the exterior edges
    std::cout << "Exterior edges:" << std::endl;
    for (const auto& edgeEntry : edgeCount) {
        if (edgeEntry.second == 1) { // Exterior edge appears only once
            std::cout << "Edge between vertices " << edgeEntry.first.first
                      << " and " << edgeEntry.first.second << std::endl;
        }
    }

    return 0;
}

int test2(){
    // Set to store unique edges and a set to store duplicate edges
    std::set<std::pair<int, int>> uniqueEdges;
    std::set<std::pair<int, int>> duplicateEdges;

    // Iterate over each triangle and its edges
    for (const auto& triangle : triangles) {
        for (int i = 0; i < 3; ++i) {
            int v1 = triangle[i];
            int v2 = triangle[(i + 1) % 3]; // Next vertex in the triangle
            std::pair<int, int> edge = createEdge(v1, v2);

            // If the edge is already in uniqueEdges, move it to duplicateEdges
            if (uniqueEdges.find(edge) != uniqueEdges.end()) {
                duplicateEdges.insert(edge);
            } else {
                uniqueEdges.insert(edge);
            }
        }
    }
}


/*

The choice between the two code implementations depends on priorities such as readability, simplicity, and performance:

std::map version:

Pros: Clear and straightforward. You directly use the map to count the number of times each edge appears, which makes it easy to identify exterior edges by checking if the count is 1.
Cons: Slightly more complex to read than the std::set version if you're unfamiliar with counting via std::map.
std::set version:

Pros: Conceptually simple, especially if you want to emphasize uniqueness. You check edges using two sets to differentiate between unique and duplicate edges.
Cons: Less efficient than the std::map version due to the need for two sets and multiple insertions/checks. The logic might also seem less direct since you need to maintain two separate sets.

*/
