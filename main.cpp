#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <SFML/Graphics.hpp>
#include <SFML/System/Clock.hpp>
#include "Line_Segment_Intersection.h";

sf::Font font;


inline float distanceBetweenVertices(sf::Vector2f pos1, sf::Vector2f pos2) {
    float& x1 = pos1.x;
    float& x2 = pos2.x;
    float& y1 = pos1.y;
    float& y2 = pos2.y;
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
}

int SGN(const float x) {
    return (x > 0) - (x < 0);
}

inline float dotProduct(const sf::Vector2f a, const sf::Vector2f b)
{
    return a.x * b.x + a.y * b.y;
}

float crossProduct(const sf::Vector2f a, const sf::Vector2f b) {
    return (a.x * b.y) - (a.y * b.x);
}

// normalizes two vectors
// see: https://stackoverflow.com/questions/10095524/normalize-a-vector
void NormalizeVector(sf::Vector2f& vec) {
    float sum = vec.x + vec.y;
    vec.x / sum;
    vec.y / sum;
}

float clamp(float val, float LB, float UB) {
    if (val < LB) {
        val = LB;
    }
    else if (val > UB) {
        val = UB;
    }
    return val;
}

namespace triangulation {

    enum WindingOrder
    {
        Clockwise, CounterClockwise, Invalid
    };

       
    float FindPolygonArea(std::vector<sf::Vector2f> vertices)
    {
        float totalArea = 0;

        for (int i = 0; i < vertices.size(); i++)
        {
            sf::Vector2f a = vertices[i];
            sf::Vector2f b = vertices[(i + 1) % vertices.size()];

            float dy = (a.y + b.y) / 2.f;
            float dx = b.x - a.x;

            float area = dy * dx;
            totalArea += area;
        }

        return abs(totalArea);
    }

    bool IsPointInTriangle(sf::Vector2f p, sf::Vector2f a, sf::Vector2f b, sf::Vector2f c)
    {
        sf::Vector2f ab = b - a;
        sf::Vector2f bc = c - b;
        sf::Vector2f ca = a - c;

        sf::Vector2f ap = p - a;
        sf::Vector2f bp = p - b;
        sf::Vector2f cp = p - c;

        float cross1 = crossProduct(ab, ap);
        float cross2 = crossProduct(bc, bp);
        float cross3 = crossProduct(ca, cp);

        if (cross1 > 0.f || cross2 > 0.f || cross3 > 0.f)
        {
            return false;
        }

        return true;
    }
    bool Triangulate(std::vector<sf::Vector2f> vertices, std::vector<int> &triangles /*the indices of the vertices for the triangles?*/)
    {

        if (vertices.size() < 3)
        {
            //errorMessage = "The vertex list must have at least 3 vertices.";
            return false;
        }

        if (vertices.size() > 1024)
        {
            //errorMessage = "The max vertex list length is 1024";
            return false;
        }

        //if (!PolygonHelper.IsSimplePolygon(vertices))
        //{
        //    errorMessage = "The vertex list does not define a simple polygon.";
        //    return false;
        //}

        //if(PolygonHelper.ContainsColinearEdges(vertices))
        //{
        //    errorMessage = "The vertex list contains colinear edges.";
        //    return false;
        //}

        //PolygonHelper.ComputePolygonArea(vertices, out float area, out WindingOrder windingOrder);

        //if(windingOrder is WindingOrder.Invalid)
        //{
        //    errorMessage = "The vertices list does not contain a valid polygon.";
        //    return false;
        //}

        //if(windingOrder is WindingOrder.CounterClockwise)
        //{
        //    Array.Reverse(vertices);
        //}

        std::vector<int> indexList;
        for (int i = 0; i < vertices.size(); i++)
        {
            indexList.push_back(i);
        }

        const int totalTriangleCount = vertices.size() - 2;
        const int totalTriangleIndexCount = totalTriangleCount * 3;

        //std::vector<int> triangles;
        triangles.resize(totalTriangleIndexCount);
        int triangleIndexCount = 0;

        while (indexList.size() > 3)
        {
            for (int i = 0; i < indexList.size(); i++)
            {
                int a = indexList[i % indexList.size()];
                int b = indexList[(i - 1) % indexList.size()];
                int c = indexList[(i + 1) % indexList.size()];


                sf::Vector2f va = vertices[a];
                sf::Vector2f vb = vertices[b];
                sf::Vector2f vc = vertices[c];

                sf::Vector2f va_to_vb = vb - va;
                sf::Vector2f va_to_vc = vc - va;

                // Is ear test vertex convex?
                if (crossProduct(va_to_vb, va_to_vc) < 0)
                {
                    continue;
                }

                bool isEar = true;

                // Does test ear contain any polygon vertices?
                for (int j = 0; j < vertices.size(); j++)
                {
                    if (j == a || j == b || j == c)
                    {
                        continue;
                    }

                    sf::Vector2f p = vertices[j];

                    if (IsPointInTriangle(p, vb, va, vc))
                    {
                        isEar = false;
                        break;
                    }
                }

                if (isEar)
                {
                    triangles[triangleIndexCount++] = b;
                    triangles[triangleIndexCount++] = a;
                    triangles[triangleIndexCount++] = c;

                    indexList.erase(indexList.begin() + i);
                    break;
                }
            }
        }

        triangles[triangleIndexCount++] = indexList[0];
        triangles[triangleIndexCount++] = indexList[1];
        triangles[triangleIndexCount++] = indexList[2];

        return true;
    }
};

namespace physics {
    #define GRAVITY .2 // 9.8 M/S

    struct SimplePolygon;
    struct ComplexPolygon;
    struct PhysicsConstraintVtoV;
    struct PhysicsConstraintCtoV;
    struct Circle;
    struct Edge;

    // a vector divided by its magnitude
    inline sf::Vector2f normalize(sf::Vector2f vec) {
        float mag = hypot(vec.x, vec.y);
        return vec / mag;
    }

    // store all sorts of shapes in here
    class physicsGroup {
        std::vector<physicsGroup> collide;
        std::vector<physicsGroup> overlap;
        std::vector<physicsGroup> ignore;

        std::vector<SimplePolygon> SimplePolygons;
        std::vector<ComplexPolygon> ComplexPolygons;
        std::vector<PhysicsConstraintVtoV> PhysicsConstraints;
        std::vector<Circle> circles;

        // add shape function
    };

    std::vector<SimplePolygon> SimplePolygons;
    std::vector<ComplexPolygon> ComplexPolygons;
    std::vector<PhysicsConstraintVtoV> PhysicsConstraintsVtoV;
    std::vector<PhysicsConstraintCtoV> PhysicsConstraintsCtoV;
    std::vector<Circle> circles;

    struct Vertex {
        sf::Vector2f Position;
        sf::Vector2f LastPos;
        sf::Vector2f Velocity = sf::Vector2f(0,0);
        float mass = 1;
        float elasticity = .1;
        void addImpulse(sf::Vector2f impulse, bool clampVelocity = false, sf::Vector2f totalVeloc = sf::Vector2f(0,0)) {
            totalVeloc.x /= 2;
            totalVeloc.y /= 2;
            if (clampVelocity) {
                impulse.x = clamp(impulse.x, -abs(totalVeloc.x), abs(totalVeloc.x));
                impulse.y = clamp(impulse.y, -abs(totalVeloc.y), abs(totalVeloc.y));
            }
            Velocity += impulse;
        }
        void invertVelocity() {
            Velocity.x *= -1;
            Velocity.y *= -1;
        }
    };

    struct Edge {
        Vertex* v1 = NULL;
        Vertex* v2 = NULL;
        
        bool Collision = true;

        float Stiffness = 1;
        float OriginalLength; // The length of the edge when it was created
        SimplePolygon* Parent = NULL; // The simple polygon that it belongs to

        sf::Vector2f getNormal() {
            sf::Vector2f temp = v1->Position - v2->Position;
            sf::Vector2f Normal = sf::Vector2f(temp.y, temp.x);
            Normal = sf::Vector2f(v1->Position.y - v2->Position.y, v2->Position.x - v1->Position.x);
            Normal = normalize(Normal);
            return Normal;
        }
        
        // returns angle in radians
        float getAngleRadians() {
            return atan2(v2->Position.y - v1->Position.y, v2->Position.x - v1->Position.x);
        }
    };

    struct SimplePolygon {
        bool anchored = false; // If it's anchored then it will not simulate physics
        bool rigidBody = true;
        float mass = 1000; // In kg

        sf::Vector2f Center;
        std::vector<Vertex> vertices;
        std::vector<Edge> edges;
        std::vector<Edge> constrainingEdges; // edges made to stop the shape from deforming

        void project_to_axis(sf::Vector2f& Axis, float& Min, float& Max) {
            float DotP = dotProduct(Axis, vertices[0].Position);

            //Set the minimum and maximum values to the projection of the first vertex
            Min = Max = DotP;

            for (int I = 1; I < vertices.size(); I++) {
                //Project the rest of the vertices onto the axis and extend
                //the interval to the left/right if necessary
                DotP = dotProduct(Axis, vertices[I].Position);

                Min = std::min(DotP, Min);
                Max = std::max(DotP, Max);
            }
        }

        // Needs to be adjusted to go off of surface area instead.
        void inline updateCenter() {
            float avgX = 0;
            float avgY = 0;
            for (int i = 0; i < vertices.size(); i++) {
                avgX += vertices[i].Position.x;
                avgY += vertices[i].Position.y;
            }
            avgX /= vertices.size();
            avgY /= vertices.size();
            Center.x = avgX;
            Center.y = avgY;
        }

        void offset(sf::Vector2f pos) {
            //todo;
        }
    };
    
    void project_to_axis(std::vector<Edge> edges, sf::Vector2f& Axis, float& Min, float& Max) {
        float DotP = dotProduct(Axis, edges[0].v1->Position);

        //Set the minimum and maximum values to the projection of the first vertex
        Min = Max = DotP;

        for (int i = 1; i < edges.size(); i++) {
            //Project the rest of the vertices onto the axis and extend
            //the interval to the left/right if necessary

            DotP = dotProduct(Axis, edges[i].v1->Position);

            Min = std::min(DotP, Min);
            Max = std::max(DotP, Max);
        }
    }

    struct ComplexPolygon {
        bool anchored = false;
        bool rigidBody = true;
        float mass;

        sf::Vector2f Center;
        sf::Vector2f Pos;
        sf::Vector2f Velocity = sf::Vector2f(0, 0);

        std::vector<Vertex> vertices;
        std::vector<Edge> Edges;
        std::vector<std::vector<Edge>> triangles;

        void project_to_axis(std::vector<Vertex*> tVerts, sf::Vector2f& Axis, float& Min, float& Max) {
            float DotP = dotProduct(Axis, tVerts.at(0)->Position);

            //Set the minimum and maximum values to the projection of the first vertex
            Min = Max = DotP;

            for (int I = 1; I < tVerts.size(); I++) {
                //Project the rest of the vertices onto the axis and extend
                //the interval to the left/right if necessary
                DotP = dotProduct(Axis, tVerts.at(I)->Position);

                Min = std::min(DotP, Min);
                Max = std::max(DotP, Max);
            }
        }

        // Needs to be adjusted to go off of surface area instead.
        sf::Vector2f updateCenter(std::vector<Vertex> triangle) {
            float avgX = 0;
            float avgY = 0;
            for (int i = 0; i < triangle.size(); i++) {
                avgX += triangle[i].Position.x;
                avgY += triangle[i].Position.y;
            }
            avgX /= triangle.size();
            avgY /= triangle.size();
            return sf::Vector2f(avgX, avgY);
        }

        void offset(sf::Vector2f pos) {
            //todo;
        }
    };

    struct PhysicsConstraintVtoV {
        float SC = 0; // stiffness constant - lower = slower movement - higher = faster movement
        float FC = 0; // friction constant (0-1)
        float mass = 1;
        float currentLength = .5;
        float targetLength = .5;
        float maxForce = 500; // the max force that can be applied by the physics constraint.
        sf::Vector2f end0Velocity;
        sf::Vector2f end1Velocity;

        Vertex* v1 = NULL;
        Vertex* v2 = NULL;
        Edge* Normal; // Edge used to calculate normal
        std::vector<float> linearRangeX = {-25,25};
        bool relative = true;
        std::vector<float> linearRangeY = {-25,25};
    };

    struct PhysicsConstraintCtoV {
        float SC = 0; // stiffness constant - lower = slower movement - higher = faster movement
        float FC = 0; // friction constant (0-1)
        float mass = 1;
        float currentLength = .5;
        float targetLength = .5;
        float maxForce = 500; // the max force that can be applied by the physics constraint.
        sf::Vector2f end0Velocity;
        sf::Vector2f end1Velocity;

        Vertex* v = NULL;
        Circle* c = NULL;
        Edge* Normal; // Edge used to calculate normal
        std::vector<float> linearRangeX = { -25,25 };
        bool relative = true;
        std::vector<float> linearRangeY = { -25,25 };
    };

    void makePhysicsConstraint(Vertex& v1, Vertex& v2, Edge &Normal, float stiffness = .5, float friction = .05, bool relative = true) {
        PhysicsConstraintVtoV PC;
        PC.targetLength = distanceBetweenVertices(v1.Position, v2.Position);
        PC.v1 = &v1;
        PC.v2 = &v2;
        PC.SC = stiffness;
        PC.FC = friction;
        if (relative) {
            PC.Normal = &Normal;
        }
        PC.relative = relative;
        PhysicsConstraintsVtoV.push_back(PC);
    }

    // Needs to be adjusted to go off of surface area instead.
    sf::Vector2f calculateCenter(std::vector<Edge> triangle) {
        float avgX = 0;
        float avgY = 0;
        for (int i = 0; i < triangle.size(); i++) {
            avgX += triangle[i].v1->Position.x;
            avgY += triangle[i].v2->Position.y;
        }
        avgX /= triangle.size();
        avgY /= triangle.size();
        return sf::Vector2f(avgX, avgY);
    }

    struct Circle {
        bool Colliding = false;


        sf::Vector2f Pos;
        sf::Vector2f Velocity = sf::Vector2f(0,0);
        float mass;
        float elasticity = 0.1; // 0 to 1
        float radius = 5; 
        void addImpulse(sf::Vector2f impulse, bool clampVelocity = false, sf::Vector2f totalVel = sf::Vector2f(0,0)) {
            totalVel.x /= 2;
            totalVel.y /= 2;
            if (clampVelocity) {
                impulse.x = clamp(impulse.x, -abs(totalVel.x), abs(totalVel.x));
                impulse.y = clamp(impulse.y, -abs(totalVel.y), abs(totalVel.y));
            }
            Velocity += (impulse);
        }
        void invertVelocity() {
            Velocity.x *= -1;
            Velocity.y *= -1;
        }
    };

    struct {
        float Depth = 0;
        sf::Vector2f Normal;

        Edge *E;
        Vertex *V;
    } CollisionInfo;


    void MakePolygon(std::vector<sf::Vector2f> Points, bool anchored_ = false) {

        SimplePolygon body_;
        SimplePolygons.push_back(body_);
        SimplePolygon &body = SimplePolygons[SimplePolygons.size() - 1];

        body.anchored = anchored_;
        // vertices
        for (int i = 0; i < Points.size(); i++) {
            Vertex V;
            V.Position = Points[i];
            V.LastPos = Points[i];
            body.vertices.push_back(V);
        }

        // edges
        // the initial edge
        Edge E;
        E.v1 = &body.vertices[0];
        E.v2 = &body.vertices[1];
        E.OriginalLength = distanceBetweenVertices(E.v1->Position, E.v2->Position);
        E.Parent = &body;
        body.edges.push_back(E);

        if (Points.size() > 2) {
            for (short int i = 1; i < Points.size() - 1; i++) {
                E.v1 = &body.vertices[i];
                E.v2 = &body.vertices[i + 1];
                E.OriginalLength = distanceBetweenVertices(E.v1->Position, E.v2->Position);
                E.Parent = &body;
                body.edges.push_back(E);
            }

            // the last edge
            E.v1 = &body.vertices[Points.size() - 1];
            E.v2 = &body.vertices[0];
            E.OriginalLength = distanceBetweenVertices(E.v1->Position, E.v2->Position);
            E.Parent = &body;
            body.edges.push_back(E);

            // the last edge
            E.v1 = &body.vertices[1];
            E.v2 = &body.vertices[3];
            E.OriginalLength = distanceBetweenVertices(E.v1->Position, E.v2->Position);
            E.Parent = &body;
            body.constrainingEdges.push_back(E);

            // the last edge
            E.v1 = &body.vertices[2];
            E.v2 = &body.vertices[0];
            E.OriginalLength = distanceBetweenVertices(E.v1->Position, E.v2->Position);
            E.Parent = &body;
            body.constrainingEdges.push_back(E);
            // the last edge
            /*E.v1 = &body.vertices[1];
            E.v2 = &body.vertices[4];
            E.OriginalLength = distanceBetweenVertices(E.v1->Position, E.v2->Position);
            E.Parent = &body;
            body.constrainingEdges.push_back(E);*/
        }
    }
    
    // set edge collision to false upon merge, but do not delete edge
    void mergeEdges(std::vector<std::vector<Edge>> &triangles) {
        for (unsigned int i = 0; i < triangles.size(); i++) {
            std::vector<Edge> &t1 = triangles[i];
            for (unsigned int j = 0; j < triangles.size(); j++) {
                std::vector<Edge>& t2 = triangles[j];
                if (i != j) {
                    for (unsigned int e1i = 0; e1i < 3; e1i++) {
                        Edge &e1 = t1[e1i];
                        for (unsigned int e2i = 0; e2i < 3; e2i++) {
                            Edge &e2 = t2[e2i];
                            if ((e1.v1 == e2.v1 && e1.v2 == e2.v2) || (e1.v1 == e2.v2 && e1.v2 == e2.v1)) {
                                if (e1.Collision) {
                                    e2.Collision = false;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void MakeComplexPolygon(std::vector<sf::Vector2f> points) {
        ComplexPolygon newPolygon;
        ComplexPolygons.push_back(newPolygon);
        ComplexPolygon &polygon = ComplexPolygons[ComplexPolygons.size() - 1];
        std::vector<Edge> &edges = polygon.Edges;
        std::vector<std::vector<Edge>> &triangles = polygon.triangles;

        std::vector<Vertex>& vertices = polygon.vertices;
        for (unsigned int i = 0; i < points.size(); i++) {
            Vertex v;
            v.mass = v.mass / vertices.size();
            v.Position = points[i];
            vertices.push_back(v);
        }
        using namespace triangulation;
        std::vector<int> indices;
       
        Triangulate(points, indices);

        for (unsigned int i = 0; i < indices.size(); i += 3) {

            std::vector<Edge> triangle (3);

            const int& a = indices[i];
            const int& b = indices[i+1];
            const int& c = indices[i+2];

            Vertex &va = vertices[a];
            Vertex &vb = vertices[b];
            Vertex &vc = vertices[c];

            Edge E;

            E.v1 = &va;
            E.v2 = &vb;
            E.OriginalLength = distanceBetweenVertices(E.v1->Position, E.v2->Position);
            triangle[0] = E;
            
            E.v1 = &vb;
            E.v2 = &vc;
            E.OriginalLength = distanceBetweenVertices(E.v1->Position, E.v2->Position);
            triangle[1] = E;

            E.v1 = &vc;
            E.v2 = &va;
            E.OriginalLength = distanceBetweenVertices(E.v1->Position, E.v2->Position);
            triangle[2] = E;

            triangles.push_back(triangle);
        }

        if (vertices.size() > 4) {
            int i = vertices.size()-1;

            std::vector<Edge> triangle(3);

            const int& a = i;
            const int& b = (i + 1)%vertices.size();
            const int& c = (i + 2)%vertices.size();

            Vertex& va = vertices[a];
            Vertex& vb = vertices[b];
            Vertex& vc = vertices[c];

            // push back first or pointer corruption will occur

            Edge E;
            E.v1 = &va;
            E.v2 = &vb;
            E.OriginalLength = distanceBetweenVertices(E.v1->Position, E.v2->Position);
            triangle[0] = E;

            E.v1 = &vb;
            E.v2 = &vc;
            E.OriginalLength = distanceBetweenVertices(E.v1->Position, E.v2->Position);
            triangle[1] = E;

            E.v1 = &vc;
            E.v2 = &va;
            E.OriginalLength = distanceBetweenVertices(E.v1->Position, E.v2->Position);
            triangle[2] = E;

            triangles.push_back(triangle);
        }

        mergeEdges(triangles);
    }

    void MakeCircle(sf::Vector2f pos, float radius) {
        Circle cir;
        cir.radius = radius;
        cir.Pos = pos;
        circles.push_back(cir);
    }
    
    // Project the vertices of each polygon onto a axis
    inline void compute_projections(const std::vector<Vertex>& bounds_a, const std::vector<Vertex>& bounds_b, const sf::Vector2f& axis_normalised, std::vector<double>& projections_a, std::vector<double>& projections_b) {
        projections_a.clear();
        projections_b.clear();

        for (size_t i = 0; i < bounds_a.size(); i++) {
            const double projection_a = dotProduct(axis_normalised, bounds_a[i].Position);


            const double projection_b = dotProduct(axis_normalised, bounds_b[i].Position);

            projections_a.push_back(projection_a);
            projections_b.push_back(projection_b);
        }
    }

    // Check if the projections of two polygons overlap
    inline bool is_overlapping(const std::vector<double>& projections_a, const std::vector<double>& projections_b) {
        const double max_projection_a = *std::max_element(projections_a.begin(), projections_a.end());
        const double min_projection_a = *std::min_element(projections_a.begin(), projections_a.end());
        const double max_projection_b = *std::max_element(projections_b.begin(), projections_b.end());
        const double min_projection_b = *std::min_element(projections_b.begin(), projections_b.end());

        // True if projection overlaps but does not necessarily mean the polygons are intersecting yet
        return !(max_projection_a < min_projection_b or max_projection_b < min_projection_a);
    }

    inline float IntervalDistance(float MinA, float MaxA, float MinB, float MaxB) {
        if (MinA < MinB)
            return MinB - MaxA;
        else
            return MinA - MaxB;
    }

    // Check if two convex polygons intersect
    inline bool separating_axis_intersect_Simple_Polygon(SimplePolygon& b1, SimplePolygon& b2, float Delta_Time) {
        b1.updateCenter();
        b2.updateCenter();

        float MinDistance = 99999.f;

        for (int i = 0; i < b1.edges.size() + b2.edges.size(); i++) {
            Edge *E;

            if (i < b1.edges.size()) {
                E = &b1.edges[i];
            }
            else {
                E = &b2.edges[i - b1.edges.size()];
            }

            sf::Vector2f Axis(E->v1->Position.y - E->v2->Position.y, E->v2->Position.x - E->v1->Position.x);
            Axis = normalize(Axis);
            
            float minA , minB, maxA , maxB; //Project both bodies onto the perpendicular axis
            b1.project_to_axis(Axis, minA, maxA);
            b2.project_to_axis(Axis, minB, maxB);

            float Distance = IntervalDistance(minA, maxA, minB, maxB);
            if (Distance > 0) {
                return false;
            }
            else {
                if (abs(Distance) < MinDistance) {
                    MinDistance = abs(Distance);
                    CollisionInfo.Normal = Axis; // invert this to switch properly colliding edges
                    CollisionInfo.E = E; //Store the edge, as it is the collision edge
                }
            }
        }

        CollisionInfo.Depth = MinDistance;

        //Ensure that the body containing the collision edge lies in
        //B2 and the one containing the collision vertex in B1
        if (CollisionInfo.E->Parent != &b2) {
            SimplePolygon Temp = b2;
            //b2 = b1;
            //b1 = Temp;
        }


        //This is needed to make sure that the collision normal is Vertexing at B1
        sf::Vector2f &BA = b1.Center;
        sf::Vector2f &BB = b2.Center;
        int Sign = SGN(dotProduct(CollisionInfo.Normal, (BA - BB)));

        //Remember that the line equation is N*( R - R0 ). We choose B2->Center
        //as R0; the normal N is given by the collision normal

        if (Sign != 1) {
            CollisionInfo.Normal = -CollisionInfo.Normal; //Revert the collision normal if it faces away from B1
        }

        float SmallestD = 99999.0f; //Initialize the smallest distance to a high value

        for (int I = 0; I < b1.vertices.size(); I++) {
            //Measure the distance of the vertex from the line using the line equation
            //float Distance = CollisionInfo.Normal * (B1->Vertices[I]->Position - B2->Center);
            float Distance = dotProduct(CollisionInfo.Normal, (b1.vertices[I].Position - BB));

            //If the measured distance is smaller than the smallest distance reported
            //so far, set the smallest distance and the collision vertex
            if (Distance < SmallestD) {
                SmallestD = Distance;
                CollisionInfo.V = &b1.vertices[I];
            }
        }
        return true;
    }

    // simple polygon to triangle SAT
    inline bool separating_axis_intersect_complex_polygon_simple_polygon(std::vector<Edge>& b1, SimplePolygon &b2, float Delta_Time) {
        float MinDistance = 99999.f;

        for (int i = 0; i < b1.size() + b2.edges.size(); i++) {
            Edge* E;

            if (i < b1.size()) {
                E = &b1[i];
            }
            else {
                E = &b2.edges[i - b1.size()];
            }

            sf::Vector2f Axis(E->v2->Position.y - E->v1->Position.y, E->v1->Position.x - E->v2->Position.x);
            Axis = normalize(Axis);

            float minA, minB, maxA, maxB; //Project both bodies onto the perpendicular axis
            project_to_axis(b1, Axis, minA, maxA);
            b2.project_to_axis(Axis, minB, maxB);

            float Distance = IntervalDistance(minA, maxA, minB, maxB);
            if (Distance > 0) {
                return false;
            }
            else {
                if (abs(Distance) < MinDistance) {
                    MinDistance = abs(Distance);
                    CollisionInfo.Normal = Axis; // invert this to switch properly colliding edges
                    CollisionInfo.E = E; //Store the edge, as it is the collision edge
                }
            }
        }

        CollisionInfo.Depth = MinDistance;

        //Ensure that the body containing the collision edge lies in
        //B2 and the one containing the collision vertex in B1
        /*if (CollisionInfo.E->Parent != b2[0].Parent) {
            std::vector<Edge> Temp = b2;
            b2 = b1;
            b1 = Temp;
        }*/


        //This is needed to make sure that the collision normal is Vertexing at B1
        sf::Vector2f BA = calculateCenter(b1);
        b2.updateCenter();
        sf::Vector2f BB = b2.Center;
        int Sign = SGN(dotProduct(CollisionInfo.Normal, (BA - BB)));

        //Remember that the line equation is N*( R - R0 ). We choose B2->Center
        //as R0; the normal N is given by the collision normal

        if (Sign != 1) {
            CollisionInfo.Normal = -CollisionInfo.Normal; //Revert the collision normal if it faces away from B1
        }

        float SmallestD = 99999.0f; //Initialize the smallest distance to a high value
        for (int I = 0; I < b1.size(); I++) {
            //Measure the distance of the vertex from the line using the line equation
            //float Distance = CollisionInfo.Normal * (B1->Vertices[I]->Position - B2->Center);
            float Distance = dotProduct(CollisionInfo.Normal, (b1[I].v1->Position - BB));

            //If the measured distance is smaller than the smallest distance reported
            //so far, set the smallest distance and the collision vertex
            if (Distance < SmallestD) {
                SmallestD = Distance;
                CollisionInfo.V = b1[I].v1;
            }
        }
        return true;
    }

    // triangle to triangle SAT
    inline bool separating_axis_intersect_triangle_to_triangle(std::vector<Edge> & b1, std::vector<Edge>& b2, float Delta_Time) {

        float MinDistance = 99999.f;

        for (int i = 0; i < b1.size() + b2.size(); i++) {
            Edge* E;

            if (i < b1.size()) {
                E = &b1[i];
            }
            else {
                E = &b2[i - b1.size()];
            }

            sf::Vector2f Axis(E->v2->Position.y - E->v1->Position.y, E->v1->Position.x - E->v2->Position.x);
            Axis = normalize(Axis);

            float minA, minB, maxA, maxB; //Project both bodies onto the perpendicular axis
            project_to_axis(b1, Axis, minA, maxA);
            project_to_axis(b2, Axis, minB, maxB);

            float Distance = IntervalDistance(minA, maxA, minB, maxB);
            if (Distance > 0) {
                return false;
            }
            else {
                if (abs(Distance) < MinDistance) {
                    MinDistance = abs(Distance);
                    CollisionInfo.Normal = Axis; // invert this to switch properly colliding edges
                    CollisionInfo.E = E; //Store the edge, as it is the collision edge
                }
            }
        }

        CollisionInfo.Depth = MinDistance;

        //Ensure that the body containing the collision edge lies in
        //B2 and the one containing the collision vertex in B1
        if (CollisionInfo.E->Parent != b2[0].Parent) {
            std::vector<Edge> Temp = b2;
            b2 = b1;
            b1 = Temp;
        }


        //This is needed to make sure that the collision normal is Vertexing at B1
        sf::Vector2f BA = calculateCenter(b1);
        sf::Vector2f BB = calculateCenter(b2);
        int Sign = SGN(dotProduct(CollisionInfo.Normal, (BA - BB)));

        //Remember that the line equation is N*( R - R0 ). We choose B2->Center
        //as R0; the normal N is given by the collision normal

        if (Sign != 1) {
            CollisionInfo.Normal = -CollisionInfo.Normal; //Revert the collision normal if it faces away from B1
        }

        float SmallestD = 99999.0f; //Initialize the smallest distance to a high value
        for (int I = 0; I < b1.size(); I++) {
            //Measure the distance of the vertex from the line using the line equation
            //float Distance = CollisionInfo.Normal * (B1->Vertices[I]->Position - B2->Center);
            float Distance = dotProduct(CollisionInfo.Normal, (b1[I].v1->Position - BB));

            //If the measured distance is smaller than the smallest distance reported
            //so far, set the smallest distance and the collision vertex
            if (Distance < SmallestD) {
                SmallestD = Distance;
                CollisionInfo.V = b1[I].v1;
            }
        }
        return true;
    }
    // polygon to circle SAT
    // https://stackoverflow.com/questions/37756135/collision-detection-separating-axis-theorem-circle-versus-polygon
    inline bool separating_axis_intersect_circle_simple_polygon(SimplePolygon &B, Circle &C, float Delta_Time) {
        
        float py = C.Pos.y;
        float px = C.Pos.x;
        bool collision = false;

        // go through each of the vertices, plus the next
        // vertex in the list
        int next = 0;
        for (int current = 0; current < B.vertices.size(); current++) {

            // get next vertex in list
            // if we've hit the end, wrap around to 0
            next = current + 1;
            if (next == B.vertices.size()) next = 0;

            // get the PVectors at our current position
            // this makes our if statement a little cleaner
            sf::Vector2f vc = B.vertices[current].Position;    // c for "current"
            sf::Vector2f vn = B.vertices[next].Position;       // n for "next"

            // compare position, flip 'collision' variable
            // back and forth
            if (((vc.y > py && vn.y < py) || (vc.y < py && vn.y > py)) &&
                (px < (vn.x - vc.x) * (py - vc.y) / (vn.y - vc.y) + vc.x)) {
                collision = !collision;
            }
        }
        // cast to edges & vertices
        float minDistance = 50000;
        Edge* CollisionEdge = NULL;
        sf::Vector2f MTV; // minimun translation vector
        for (unsigned int i = 0; i < B.edges.size(); i++) {
            Edge& E = B.edges[i];
            sf::Vector2f Normal = sf::Vector2f(E.v1->Position.y - E.v2->Position.y, E.v2->Position.x - E.v1->Position.x);
            Normal.x *= -50000;
            Normal.y *= -50000;
            sf::Vector2f p1 = C.Pos;
            sf::Vector2f p2 = C.Pos + Normal;
            sf::Vector2f pInt = CalcIntersection(p1, p2, E.v1->Position, E.v2->Position);
            float dist = distanceBetweenVertices(pInt, C.Pos);
            if (dist < C.radius) {
                if (dist < minDistance) {
                    minDistance = dist;
                    CollisionEdge = &E;
                    MTV = (pInt - C.Pos);
                }
            }
        }
        
        if (minDistance != 50000) {
            const float weightV2 = 1 / (distanceBetweenVertices(C.Pos, CollisionEdge->v2->Position) / distanceBetweenVertices(C.Pos, CollisionEdge->v1->Position));
            const float weightV1 = 1 / (distanceBetweenVertices(C.Pos, CollisionEdge->v1->Position) / distanceBetweenVertices(C.Pos, CollisionEdge->v2->Position));
            float two = 2;
            float four = 4;
            sf::Vector2f totalVel1 = sf::Vector2f(abs(C.Velocity.x) + abs(CollisionEdge->v1->Velocity.x), abs(C.Velocity.y) + abs(CollisionEdge->v1->Velocity.y));
            sf::Vector2f totalVel2 = sf::Vector2f(abs(C.Velocity.x) + abs(CollisionEdge->v2->Velocity.x), abs(C.Velocity.y) + abs(CollisionEdge->v2->Velocity.y));

            C.addImpulse(-(MTV / C.radius) / two, true, (totalVel1 + totalVel2) / two);

            CollisionEdge->v1->addImpulse(sf::Vector2f(MTV.x / weightV1, MTV.y / weightV1) / C.radius / two, true, totalVel1 * weightV1);
            CollisionEdge->v2->addImpulse(sf::Vector2f(MTV.x / weightV2, MTV.y / weightV2) / C.radius / two, true, totalVel2 * weightV2);
        
            return true;
        }

        for (unsigned int i = 0; i < B.vertices.size(); i++) {
            Vertex& V = B.vertices[i];
            sf::Vector2f diff = C.Pos - V.Position;
            if (hypot(diff.x, diff.y) < C.radius) {
                if (C.Colliding == false) {
                    float t = C.radius;
                    sf::Vector2f combinedVel = sf::Vector2f(abs(C.Velocity.x) + abs(V.Velocity.x), abs(C.Velocity.y) + abs(V.Velocity.y));
                    C.Pos += diff / t;
                    C.addImpulse(diff / t, true, combinedVel);
                    V.addImpulse(-diff / t, true, combinedVel);
                }
            }
        }
        return false;
    }

    inline bool separating_axis_intersect_circle_complex_polygon(ComplexPolygon &B, std::vector<Edge>& T, Circle& C, float Delta_Time) {
        float py = C.Pos.y;
        float px = C.Pos.x;
        bool collision = false;

        // go through each of the vertices, plus the next
        // vertex in the list
        int next = 0;
        for (int current = 0; current < B.vertices.size(); current++) {

            // get next vertex in list
            // if we've hit the end, wrap around to 0
            next = current + 1;
            if (next == B.vertices.size()) next = 0;

            // get the PVectors at our current position
            // this makes our if statement a little cleaner
            sf::Vector2f vc = B.vertices[current].Position;    // c for "current"
            sf::Vector2f vn = B.vertices[next].Position;       // n for "next"

            // compare position, flip 'collision' variable
            // back and forth
            if (((vc.y > py && vn.y < py) || (vc.y < py && vn.y > py)) &&
                (px < (vn.x - vc.x) * (py - vc.y) / (vn.y - vc.y) + vc.x)) {
                collision = !collision;
            }
        }

        // cast to edges & vertices
        float minDistance = 50000;
        Edge* CollisionEdge = NULL;
        sf::Vector2f MTV; // minimun translation vector
        for (unsigned int i = 0; i < T.size(); i++) {
            Edge* E = &T[i];
            sf::Vector2f Normal = sf::Vector2f(E->v1->Position.y - E->v2->Position.y, E->v2->Position.x - E->v1->Position.x);
            Normal.x *= -50000;
            Normal.y *= -50000;
            sf::Vector2f p1 = C.Pos;
            sf::Vector2f p2 = C.Pos + Normal;
            sf::Vector2f pInt = CalcIntersection(p1, p2, E->v1->Position, E->v2->Position);
            float dist = distanceBetweenVertices(pInt, C.Pos);
            if (dist < C.radius) {
                if (dist < minDistance) {
                    minDistance = dist;
                    CollisionEdge = E;
                    MTV = (pInt - C.Pos);
                }
            }
        }

        if (minDistance != 50000) {
            const float weightV2 = 1 / (distanceBetweenVertices(C.Pos, CollisionEdge->v2->Position) / distanceBetweenVertices(C.Pos, CollisionEdge->v1->Position));
            const float weightV1 = 1 / (distanceBetweenVertices(C.Pos, CollisionEdge->v1->Position) / distanceBetweenVertices(C.Pos, CollisionEdge->v2->Position));
            float two = 2;
            float four = 4;
            sf::Vector2f totalVel1 = sf::Vector2f(abs(C.Velocity.x) + abs(CollisionEdge->v1->Velocity.x), abs(C.Velocity.y) + abs(CollisionEdge->v1->Velocity.y));
            sf::Vector2f totalVel2 = sf::Vector2f(abs(C.Velocity.x) + abs(CollisionEdge->v2->Velocity.x), abs(C.Velocity.y) + abs(CollisionEdge->v2->Velocity.y));

            C.addImpulse(-(MTV / C.radius) / two, false, (totalVel1 + totalVel2) / two);
            CollisionEdge->v1->addImpulse(sf::Vector2f(MTV.x / weightV1, MTV.y / weightV1) / C.radius / two, false, totalVel1 * weightV1);
            CollisionEdge->v2->addImpulse(sf::Vector2f(MTV.x / weightV2, MTV.y / weightV2) / C.radius / two, false, totalVel2 * weightV2);
            return true;
        }

        for (unsigned int i = 0; i < B.vertices.size(); i++) {
            Vertex& V = B.vertices[i];
            sf::Vector2f diff = C.Pos - V.Position;
            if (hypot(diff.x, diff.y) < C.radius) {
                float t = C.radius;
                sf::Vector2f combinedVel = sf::Vector2f(abs(C.Velocity.x) + abs(V.Velocity.x), abs(C.Velocity.y) + abs(V.Velocity.y));
                C.Pos += diff / t;
                C.addImpulse(diff / t, true, combinedVel);
                V.addImpulse(-diff / t, true, combinedVel);
            }
        }
        return false;

    }

    inline void collisionResponse(float Delta_Time) {

        sf::Vector2f CollisionVector = CollisionInfo.Normal * CollisionInfo.Depth;
        Vertex* V1 = CollisionInfo.E->v1;
        Vertex* V2 = CollisionInfo.E->v2;

        float T;
        if (V1 == NULL || V2 == NULL) {
            return;
        }

        if (abs(V1->Position.x - V2->Position.x) > abs(V1->Position.y - V2->Position.y))
            T = (CollisionInfo.V->Position.x - CollisionVector.x - V1->Position.x) / (V2->Position.x - V1->Position.x);
        else
            T = (CollisionInfo.V->Position.y - CollisionVector.x - V1->Position.x) / (V2->Position.y - V1->Position.y);


        float Lambda = 1.0f / (T * T + (1 - T) * (1 - T));

        sf::Vector2f vel1 = V1->Velocity;
        sf::Vector2f vel2 = V2->Velocity;
        vel1.x = abs(vel1.x);
        vel1.y = abs(vel1.y);
        vel2.x = abs(vel2.x);
        vel2.y = abs(vel2.y);
        //sf::Vector2f V_Veloc = CollisionInfo.V->Velocity;
        //V_Veloc.x = abs(V_Veloc.x);
        //V_Veloc.y = abs(V_Veloc.y);
        //const float weightV2 = 1 / (distanceBetweenVertices(CollisionInfo.V->Position, V2->Position) / distanceBetweenVertices(CollisionInfo.V->Position, V1->Position));
        //const float weightV1 = 1 / (distanceBetweenVertices(CollisionInfo.V->Position, V2->Position) / distanceBetweenVertices(CollisionInfo.V->Position, V2->Position));
        V2->Position += CollisionVector * (1 - T) * 0.5f * Lambda;
        V1->Position += CollisionVector * T * 0.5f * Lambda;
        //V2->addImpulse(CollisionVector * (1 - T) * 0.5f * Lambda, false, vel1);
        //V1->addImpulse(CollisionVector * T * 0.5f * Lambda, false, vel2);
        
        float oneHalf = .5;
        float nineTenths = .9;

        sf::Vector2f addVel;
        addVel = (CollisionInfo.V->Velocity - ((V1->Velocity + V2->Velocity) * oneHalf)) * oneHalf;
        V1->addImpulse(addVel * nineTenths);
        V2->addImpulse(addVel * nineTenths);
        //CollisionInfo.V->addImpulse(V2->Velocity);

        CollisionInfo.V->Position += sf::Vector2f(CollisionVector.x * .5, CollisionVector.y * .5);
        //CollisionInfo.V->addImpulse(sf::Vector2f(CollisionVector.x,CollisionVector.y), false, vel1 + vel2);
        //(CollisionInfo.V->Velocity - ((V1->Velocity + V2->Velocity) * oneHalf))
        addVel = (((V1->Velocity + V2->Velocity) * oneHalf) - CollisionInfo.V->Velocity);
        CollisionInfo.V->addImpulse(addVel);



        // simulates friction
        // FrictionVector = ?*m*g;
        // g = gravitational constant (9.8 on earth)
        // .01 is friction coefficient
        float Friction = 0 * 2 * GRAVITY;
        sf::Vector2f FrictionVector(V1->Velocity.x * Friction, V1->Velocity.y * Friction);
        //V1->Position += FrictionVector;
        FrictionVector.x = V2->Velocity.x * Friction;
        FrictionVector.y = V2->Velocity.y * Friction;
        //V2->Position += FrictionVector
        
    }

    void updateEdges(SimplePolygon &body, bool Delta_Time) {
        for (unsigned int i = 0; i < body.edges.size() + body.constrainingEdges.size(); i++) {

            Edge* E;
            if (i < body.edges.size()) {
                E = &body.edges[i];
            } else {
                E = &body.constrainingEdges[i-body.edges.size()];
            }

            if (E->Collision) {
                float dx = E->v2->Position.x - E->v1->Position.x;
                float dy = E->v2->Position.y - E->v1->Position.y;
                float dist = distanceBetweenVertices(E->v1->Position, E->v2->Position);

                if (dist != 0) { // prevents divide by zero error
                    float diff = (E->OriginalLength - dist) / dist;

                    // gets offset of the Vertexs
                    float offsetx = dx * diff * 0.5;
                    float offsety = dy * diff * 0.5;

                    // calculate "mass"
                    float m1 = 2 + 1;
                    float m2 = 2 / m1;
                    m1 = 1 / m1;
                    m1 = 1;
                    m2 = 1;

                    E->v1->Position.x -= offsetx * m1;
                    E->v1->Position.y -= offsety * m1;
                    E->v2->Position.x += offsetx * m2;
                    E->v2->Position.y += offsety * m2;

                    E->v1->addImpulse(sf::Vector2f(-(offsetx * m1), -(offsety * m1)));
                    E->v2->addImpulse(sf::Vector2f(offsetx * m2, offsety * m2));
                }
            }
        }
    }

    void updateEdges(ComplexPolygon& body) {
           for (unsigned int t = 0; t < body.triangles.size(); t++) {
            std::vector<Edge> &triangle = body.triangles[t];
            for (unsigned int i = 0; i < triangle.size(); i++) {

                Edge* E;

                E = &triangle[i];

                float dx = E->v2->Position.x - E->v1->Position.x;
                float dy = E->v2->Position.y - E->v1->Position.y;
                float dist = distanceBetweenVertices(E->v1->Position, E->v2->Position);

                if (dist != 0) { // prevents divide by zero error
                    float diff = (E->OriginalLength - dist) / dist;

                    // gets offset of the Vertexs
                    float offsetx = dx * diff * 0.5;
                    float offsety = dy * diff * 0.5;

                    // calculate "mass"
                    float m1 = 2 + 1;
                    float m2 = 2 / m1;
                    m1 = 1 / m1;
                    m1 = 1;
                    m2 = 1;

                    E->v1->Position.x -= offsetx * m1;
                    E->v1->Position.y -= offsety * m1;
                    E->v2->Position.x += offsetx * m2;
                    E->v2->Position.y += offsety * m2;

                    E->v1->addImpulse(sf::Vector2f(-(offsetx * m1), -(offsety * m1)));
                    E->v2->addImpulse(sf::Vector2f(offsetx * m2, offsety * m2));
                }
            }
        }
    }

    inline void IterateCollisions(float Delta_Time) {

        // recenters shape. Does not currently account for volume. Center is also recalculated elsewhere in program. I need to fix this.
        for (short int i = 0; i < SimplePolygons.size(); i++) {
            SimplePolygon& body = SimplePolygons[i];
            // recalculates body center
            float totalX = 0;
            float totalY = 0;
            for (short int j = 0; j < body.vertices.size(); j++) {
                totalX += body.vertices[j].Position.x;
                totalY += body.vertices[j].Position.y;
            }
            totalX /= body.vertices.size();
            totalY /= body.vertices.size();
            body.Center = sf::Vector2f(totalX, totalY);
        }

        // simple polygons
        for (int i = 0; i < SimplePolygons.size(); i++) {
            SimplePolygon& B1 = SimplePolygons[i];
            for (int j = 0; j < SimplePolygons.size(); j++) {
                SimplePolygon& B2 = SimplePolygons[j];
                if (i != j && !SimplePolygons[j].anchored) {
                    if (separating_axis_intersect_Simple_Polygon(B1, B2, Delta_Time) && CollisionInfo.V != NULL) {
                        collisionResponse(Delta_Time);
                    }
                }
            }
            // polygon to circle collision
            for (unsigned int c = 0; c < circles.size(); c++) {
                if (separating_axis_intersect_circle_simple_polygon(B1, circles[c], Delta_Time)) {
                    circles[c].Colliding = true;
                }
                else {
                    circles[c].Colliding = false;
                }
            }
            updateEdges(B1, Delta_Time); // updates edges should be called inside the iterate collisions funciton
        }

        // complex polygons
        for (unsigned int i = 0; i < ComplexPolygons.size(); i++) {
            ComplexPolygon& B1 = ComplexPolygons[i];
            for (unsigned int j = 0; j < ComplexPolygons.size(); j++) {
                ComplexPolygon& B2 = ComplexPolygons[j];
                if (i != j && !B2.anchored) {
                    for (unsigned int a = 0; a < B1.triangles.size(); a++) {
                        std::vector<Edge> t1 = B1.triangles[a];
                        for (unsigned int b = 0; b < B2.triangles.size(); b++) {
                            std::vector<Edge> t2 = B2.triangles[b];
                            if (separating_axis_intersect_triangle_to_triangle(t1, t2, Delta_Time)) {
                                collisionResponse(Delta_Time);
                            }
                        }
                    }
                }
            }
            // complex polygon to circle collision
            for (unsigned int c = 0; c < circles.size(); c++) {
                for (unsigned int a = 0; a < B1.triangles.size(); a++) {
                    std::vector<Edge>& T = B1.triangles[a];
                    if (separating_axis_intersect_circle_complex_polygon(B1, T, circles[c], Delta_Time)) {

                    }
                }
            }
            // complex polygon to simple polygon collision?
            updateEdges(B1);
            for (unsigned int t = 0; t < B1.triangles.size(); t++) {
                for (unsigned int b = 0; b < SimplePolygons.size(); b++) {

                    std::vector<Edge>& T = B1.triangles[t];
                    if (separating_axis_intersect_complex_polygon_simple_polygon(T, SimplePolygons[b], Delta_Time)) {
                        collisionResponse(Delta_Time);
                    }
                }
            }
        }

        // circle to circle collision
        for (unsigned int i = 0; i < circles.size(); i++) {
            Circle &C1 = circles[i];
            for (unsigned int j = 0; j < circles.size(); j++) {
                Circle& C2 = circles[j];
                if (i != j) {
                    sf::Vector2f diff = C1.Pos - C2.Pos;
                    if (hypot(diff.x, diff.y) < C1.radius + C2.radius) {
                        float t = (C1.radius + C2.radius);
                        C1.Pos += diff / t;
                        C2.Pos -= diff / t;

                        sf::Vector2f Vel1 = C1.Velocity;
                        sf::Vector2f Vel2 = C2.Velocity;
                        Vel1.x = abs(Vel1.x);
                        Vel1.y = abs(Vel1.y);
                        Vel2.x = abs(Vel2.x);
                        Vel2.y = abs(Vel2.y);
                        C1.addImpulse(sf::Vector2f(diff.x, diff.y), true, sf::Vector2f(Vel1.x + Vel2.x, Vel1.y+Vel2.y));
                        C2.addImpulse(sf::Vector2f(-diff.x, -diff.y),true, sf::Vector2f(Vel1.x + Vel2.x, Vel1.y + Vel2.y));
                    }
                }
            }
        }
    }

    // applies velocity and gravity
    inline void UpdateForces(float Delta_Time) {
        // moves shapes by velocity
        // simple polygons
        for (unsigned int i = 0; i < SimplePolygons.size(); i++) {
            SimplePolygon &body = SimplePolygons[i];
            if (!body.anchored) {
                float gravForce = GRAVITY;
                float gravAccel = gravForce / body.mass;
                // verts
                for (int i = 0; i < body.vertices.size(); i++) {
                    Vertex &V = body.vertices[i];
                    V.Position += V.Velocity * Delta_Time;
                }
            }
        }

        // complex polygons
        for (unsigned int i = 0; i < ComplexPolygons.size(); i++) {
            ComplexPolygon &poly = ComplexPolygons[i];
            for (unsigned int v = 0; v < poly.vertices.size(); v++) {
                Vertex &V = poly.vertices[v];
                V.Position += V.Velocity * Delta_Time;
            }
        }

        // circles
        for (unsigned int i = 0; i < circles.size(); i++) {
            Circle& c = circles[i];
            c.Pos += c.Velocity * Delta_Time;
        }
    }

    // angle must be in radians, not degrees
    sf::Vector2f rotateVector(sf::Vector2f v, float angle) {
        float s = sin(angle);
        float c = cos(angle);

        // rotate point
        float xnew = v.x * c - v.y * s;
        float ynew = v.x * s + v.y * c;

        return sf::Vector2f(xnew, ynew);
    }

    void updatePhysicsConstraints(float dt) {
        for (unsigned int i = 0; i < PhysicsConstraintsVtoV.size(); ++i) {
            PhysicsConstraintVtoV &PC = PhysicsConstraintsVtoV[i];
            sf::Vector2f& end0 = PC.v1->Position;
            sf::Vector2f& end1 = PC.v2->Position;
            float diffx = end1.x - end0.x;
            float diffy = end1.y - end0.y;
            sf::Vector2f diff = end1 - end0;

            float length = distanceBetweenVertices(end0, end1);
            float displacement = length - PC.targetLength;

            sf::Vector2f multiplier = sf::Vector2f(1, 1);

            sf::Vector2f springN;
            if (length < .0001) { // prevents divide by zero error
                springN.x = diff.x / .0001;
                springN.y = diff.y / .0001;
            }
            else {
                springN.x = diff.x / length;
                springN.y = diff.y / length;
            }

            sf::Vector2f restoreForce = springN * (displacement * PC.SC);
            restoreForce.x *= abs(multiplier.x);
            restoreForce.y *= abs(multiplier.y);

            sf::Vector2f end0Force;
            sf::Vector2f end1Force;

            end0Force.x = restoreForce.x * .4 / PC.mass;
            end0Force.y = restoreForce.y * .4 / PC.mass;
            end1Force.x = -restoreForce.x * .4 / PC.mass;
            end1Force.y = -restoreForce.y * .4 / PC.mass;

            float fric = (1.0 - PC.FC);
            PC.end0Velocity += end0Force * dt;
            PC.end1Velocity += end1Force * dt;
            PC.end0Velocity.x *= fric;
            PC.end0Velocity.y *= fric;
            PC.end1Velocity.x *= fric;
            PC.end1Velocity.y *= fric;
            PC.v1->addImpulse(end0Force);
            PC.v2->addImpulse(end1Force);

            float minX = PC.linearRangeX[0];
            float maxX = PC.linearRangeX[1];

            float minY = PC.linearRangeY[0];
            float maxY = PC.linearRangeY[1];

            // limits constraint
            if (diffx < minX) {
                PC.v1->Position += sf::Vector2f((diffx - minX) * multiplier.x, 0);
                PC.v2->Position -= sf::Vector2f((diffx - minX) * multiplier.x, 0);
            }
            else if (diffx > maxX) {
                PC.v1->Position += sf::Vector2f((diffx - maxX) * multiplier.x, 0);
                PC.v2->Position -= sf::Vector2f((diffx - maxX) * multiplier.x, 0);
            }
            if (diffy < minY) {
                PC.v1->Position += sf::Vector2f(0, (diffy - minY) * multiplier.y);
                PC.v2->Position -= sf::Vector2f(0, (diffy - minY) * multiplier.y);
            }
            else if (diffy > maxY) {
                PC.v1->Position += sf::Vector2f(0, (diffy - maxY) * multiplier.y);
                PC.v2->Position -= sf::Vector2f(0, (diffy - maxY) * multiplier.y);
            }
            
            if (PC.relative) {
                const int two = 2;
                const sf::Vector2f normal = PC.Normal->getNormal();
                sf::Vector2f begin = (PC.v1->Position + PC.v2->Position); begin.x /= 2, begin.y /= 2;
                sf::Vector2f end = (PC.v1->Position + PC.v2->Position); end.x /= 2, end.y /= 2;
                end += normal * length;

                float restoreForce = (displacement * PC.SC);
                // limits constraint
                sf::Vector2f translation = (distanceBetweenVertices(PC.v1->Position, PC.v2->Position) * normal);
                translation.x = clamp(translation.x, -25, 25);
                translation.y = clamp(translation.y, -25, 25);
                PC.v2->Position = PC.v1->Position + translation;
            }
        }




        for (unsigned int i = 0; i < PhysicsConstraintsCtoV.size(); ++i) {
            PhysicsConstraintCtoV &PC = PhysicsConstraintsCtoV[i];
            sf::Vector2f& end0 = PC.v->Position;
            sf::Vector2f& end1 = PC.c->Pos;
            float diffx = end1.x - end0.x;
            float diffy = end1.y - end0.y;
            sf::Vector2f diff = end1 - end0;

            float length = distanceBetweenVertices(end0, end1);
            float displacement = length - PC.targetLength;

            sf::Vector2f multiplier = sf::Vector2f(1, 1);

            sf::Vector2f springN;
            if (length < .0001) { // prevents divide by zero error
                springN.x = diff.x / .0001;
                springN.y = diff.y / .0001;
            }
            else {
                springN.x = diff.x / length;
                springN.y = diff.y / length;
            }

            sf::Vector2f restoreForce = springN * (displacement * PC.SC);
            restoreForce.x *= abs(multiplier.x);
            restoreForce.y *= abs(multiplier.y);

            sf::Vector2f end0Force;
            sf::Vector2f end1Force;

            end0Force.x = restoreForce.x * .4 / PC.mass;
            end0Force.y = restoreForce.y * .4 / PC.mass;
            end1Force.x = -restoreForce.x * .4 / PC.mass;
            end1Force.y = -restoreForce.y * .4 / PC.mass;

            float fric = (1.0 - PC.FC);
            PC.end0Velocity += end0Force * dt;
            PC.end1Velocity += end1Force * dt;
            PC.end0Velocity.x *= fric;
            PC.end0Velocity.y *= fric;
            PC.end1Velocity.x *= fric;
            PC.end1Velocity.y *= fric;
            PC.v->addImpulse(end0Force);
            PC.c->addImpulse(end1Force);

            float minX = PC.linearRangeX[0];
            float maxX = PC.linearRangeX[1];

            float minY = PC.linearRangeY[0];
            float maxY = PC.linearRangeY[1];

            // limits constraint
            if (diffx < minX) {
                PC.v->Position += sf::Vector2f((diffx - minX) * multiplier.x, 0);
                PC.c->Pos -= sf::Vector2f((diffx - minX) * multiplier.x, 0);
            }
            else if (diffx > maxX) {
                PC.v->Position += sf::Vector2f((diffx - maxX) * multiplier.x, 0);
                PC.c->Pos -= sf::Vector2f((diffx - maxX) * multiplier.x, 0);
            }
            if (diffy < minY) {
                PC.v->Position += sf::Vector2f(0, (diffy - minY) * multiplier.y);
                PC.c->Pos -= sf::Vector2f(0, (diffy - minY) * multiplier.y);
            }
            else if (diffy > maxY) {
                PC.v->Position += sf::Vector2f(0, (diffy - maxY) * multiplier.y);
                PC.c->Pos -= sf::Vector2f(0, (diffy - maxY) * multiplier.y);
            }
            
            if (PC.relative) {
                const int two = 2;
                const sf::Vector2f normal = PC.Normal->getNormal();
                sf::Vector2f begin = (PC.v->Position + PC.v->Position); begin.x /= 2, begin.y /= 2;
                sf::Vector2f end = (PC.c->Pos + PC.c->Pos); end.x /= 2, end.y /= 2;
                end += normal * length;

                float restoreForce = (displacement * PC.SC);
                // limits constraint
                sf::Vector2f translation = (distanceBetweenVertices(PC.v->Position, PC.c->Pos) * normal);
                translation.x = clamp(translation.x, -25, 25);
                translation.y = clamp(translation.y, -25, 25);
                PC.c->Pos = PC.v->Position + translation;
            }
        }
    }

    inline void update(float Delta_Time) {

        // call update edges multiple times to maintain shape
        // note that when an edge updates it may interfere with the length of another

        UpdateForces(Delta_Time); // applies velocity & gravity

        updatePhysicsConstraints(Delta_Time);

        for (unsigned int i = 0; i < 50; ++i) {
            IterateCollisions(Delta_Time);
        }

        // updates the velocity of shape vertices
        for (unsigned int i = 0; i < SimplePolygons.size(); i++) {
            SimplePolygon& P = SimplePolygons[i];
            for (unsigned int j = 0; j < P.vertices.size(); j++) {
                Vertex& V = P.vertices[j];
                //V.Velocity = (V.Position - V.LastPos) / Delta_Time;
                V.LastPos = V.Position;
            }
        }
    }

    // drags vertex when user is holding down left click on vertex
    bool debuffLeft = false;
    bool usingCircle = false;
    Vertex* closestV = NULL;
    Circle* closestC = NULL;
    int closestDistV = 500;
    int closestDistC = 500;

    void dragShapes(sf::Vector2f mousePos) {
        if (!debuffLeft) {
            closestDistV = 500;
            for (unsigned int b = 0; b < ComplexPolygons.size(); b++) {
                ComplexPolygon& body = ComplexPolygons[b];
                for (int i = 0; i < body.vertices.size(); i++) {
                    if (distanceBetweenVertices(mousePos, body.vertices[i].Position) < closestDistV) {
                        closestDistV = distanceBetweenVertices(mousePos, body.vertices[i].Position);
                        closestV = &body.vertices[i];
                    }
                }
            }

            for (unsigned int b = 0; b < SimplePolygons.size(); b++) {
                SimplePolygon& body = SimplePolygons[b];
                for (int i = 0; i < body.vertices.size(); i++) {
                    if (distanceBetweenVertices(mousePos, body.vertices[i].Position) < closestDistV) {
                        closestDistV = distanceBetweenVertices(mousePos, body.vertices[i].Position);
                        closestV = &body.vertices[i];
                    }
                }
            }
            if (closestDistV < 25) {
                usingCircle = false;
            }
            closestDistC = 500;
            for (unsigned int c = 0; c < circles.size(); c++) {
                if (distanceBetweenVertices(mousePos, circles[c].Pos) < closestDistC) {
                    closestDistC = distanceBetweenVertices(mousePos, circles[c].Pos);
                    closestC = &circles[c];
                }
            }
        }

        if (closestDistV < 25) {
            const float five = 5;
            closestV->addImpulse((mousePos-closestV->Position)/five);
            //closestV->Position = mousePos;
            usingCircle = false;
        }
        if (closestC != NULL) {
            if (closestDistC < closestC->radius) {
                const float five = 5;
                closestC->addImpulse((mousePos-closestC->Pos)/five);
                usingCircle = true;
            }
        }
        debuffLeft = true;
    }

    void showVertexIndex(sf::Vector2f mousePos, sf::RenderWindow &window) {
        int closestDist = 500;
        int closestIndex;
        Vertex* closest = NULL;
        for (unsigned int b = 0; b < SimplePolygons.size(); b++) {
            SimplePolygon& body = SimplePolygons[b];
            for (int i = 0; i < body.vertices.size(); i++) {
                if (distanceBetweenVertices(mousePos, body.vertices[i].Position) < closestDist) {
                    closestIndex = i;
                    closestDist = distanceBetweenVertices(mousePos, body.vertices[i].Position);
                    closest = &body.vertices[i];
                }
            }
        }
        if (closestDist < 50) {
            sf::Text text;
            text.setFont(font);
            text.setCharacterSize(25);
            text.setPosition(sf::Vector2f(10, 10));
            text.setFillColor(sf::Color::White);
            text.setString(std::to_string(closestIndex));
            window.draw(text);
        }
    }

    void contrainShapes(std::vector<SimplePolygon>& polygons, sf::RenderWindow& window) {
        int winX = window.getSize().x;
        int winY = window.getSize().y;
        
        for (unsigned int b = 0; b < polygons.size(); b++) {
            SimplePolygon& body = polygons[b];
            for (int i = 0; i < body.vertices.size(); i++) {
                bool hit = false;
                Vertex &V = body.vertices[i];
                sf::Vector2f& pos = V.Position;
                if (pos.x < 0) { 
                    pos.x = 0; V.Velocity.x *= -1;
                }
                if (pos.y < 0) { 
                    pos.y = 0; V.Velocity.y *= -1;;
                }
                if (pos.x > winX) { 
                    pos.x = winX; V.Velocity.x *= -1;
                }
                if (pos.y > winY) { 
                    pos.y = winY; V.Velocity.y *= -1;
                }
            }
        }

        for (unsigned int b = 0; b < ComplexPolygons.size(); b++) {
            ComplexPolygon& body = ComplexPolygons[b];
            for (int i = 0; i < body.vertices.size(); i++) {
                Vertex& V = body.vertices[i];
                sf::Vector2f& pos = V.Position;
                if (pos.x < 0) {
                    pos.x = 0; V.Velocity.x *= -1;
                }
                if (pos.y < 0) {
                    pos.y = 0; V.Velocity.y *= -1;
                }
                if (pos.x > winX) {
                    pos.x = winX; V.Velocity.x *= -1;
                }
                if (pos.y > winY) {
                    pos.y = winY; V.Velocity.y *= -1;
                }
            }
        }

        for (unsigned int c = 0; c < circles.size(); c++) {
            sf::Vector2f &Pos = circles[c].Pos;
            float& r = circles[c].radius;
            Circle& circle = circles[c];
            bool hit = false;
            if (Pos.x - r < 0) {
                Pos += sf::Vector2f(0 - (Pos.x - r), 0);
                circle.Velocity.x *= -1;
            }
            if (Pos.x + r > winX) {
                Pos += sf::Vector2f(winX - (Pos.x + r), 0);
                circle.Velocity.x *= -1;
            }
            if (Pos.y - r < 0) {
                Pos += sf::Vector2f(0, 0 - (Pos.y - r));
                circle.Velocity.y *= -1;
            }
            if (Pos.y + r > winY) {
                Pos += sf::Vector2f(0, winY - (Pos.y + r));
                circle.Velocity.y *= -1;
            }
        }
    }
}

void renderVertices(sf::RenderWindow& window) {
    using namespace std;
    using namespace physics;
    // polygons


    for (int i = 0; i < physics::SimplePolygons.size(); i++) {
        physics::SimplePolygon &body = physics::SimplePolygons[i];
        for (int j = 0; j < body.vertices.size(); j++) {
            physics::Vertex V = body.vertices[j];
            sf::CircleShape circle(4);
            circle.setPosition(V.Position);
            circle.setOrigin(circle.getGlobalBounds().width / 2, circle.getGlobalBounds().height / 2);
            circle.setFillColor(sf::Color::Red);
            window.draw(circle);
        }

        for (int j = 0; j < body.edges.size(); j++) {
            physics::Vertex V1 = *body.edges[j].v1;
            sf::CircleShape circle(3);
            circle.setPosition(V1.Position);
            circle.setOrigin(circle.getGlobalBounds().width / 2, circle.getGlobalBounds().height / 2);
            circle.setFillColor(sf::Color::Green);
            window.draw(circle);

            physics::Vertex V2 = *body.edges[j].v2;
            circle.setPosition(V2.Position);
            circle.setFillColor(sf::Color::Blue);
            window.draw(circle);

            sf::RectangleShape rect(sf::Vector2f(2, distanceBetweenVertices(V1.Position, V2.Position)));
            rect.setOrigin(rect.getGlobalBounds().width / 2, rect.getGlobalBounds().height / 2);
            float avgX = (V1.Position.x + V2.Position.x) / 2;
            float avgY = (V1.Position.y + V2.Position.y) / 2;
            float rotation = atan2(V1.Position.y - V2.Position.y, V1.Position.x - V2.Position.x);
            // converts radians to degrees
            rotation = rotation * 180 / 3.14159;
            rect.setPosition(avgX, avgY);
            rect.setRotation(rotation + 90);
            rect.setFillColor(sf::Color::Green);
            window.draw(rect);

        }

        for (int j = 0; j < body.constrainingEdges.size(); j++) {
            physics::Vertex V1 = *body.constrainingEdges[j].v1;

            physics::Vertex V2 = *body.constrainingEdges[j].v2;
         

            sf::RectangleShape rect(sf::Vector2f(2, distanceBetweenVertices(V1.Position, V2.Position)));
            rect.setOrigin(rect.getGlobalBounds().width / 2, rect.getGlobalBounds().height / 2);
            float avgX = (V1.Position.x + V2.Position.x) / 2;
            float avgY = (V1.Position.y + V2.Position.y) / 2;
            float rotation = atan2(V1.Position.y - V2.Position.y, V1.Position.x - V2.Position.x);
            // converts radians to degrees
            rotation = rotation * 180 / 3.14159;
            rect.setPosition(avgX, avgY);
            rect.setRotation(rotation + 90);
            rect.setFillColor(sf::Color::Blue);
            window.draw(rect);
        }

        // draws a line showing the collision vector
        physics::Vertex V1;
        physics::Vertex V2;
        V1.Position = sf::Vector2f(0, 0);
        V2.Position = physics::CollisionInfo.Normal*physics::CollisionInfo.Depth;
        float temp = V2.Position.x;
        V2.Position.x = -V2.Position.y;
        V2.Position.y = temp;
        V2.Position.x *= 5;
        V2.Position.y *= 5;
        sf::RectangleShape rect(sf::Vector2f(2, distanceBetweenVertices(V1.Position, V2.Position)));
        rect.setOrigin(rect.getGlobalBounds().width / 2, rect.getGlobalBounds().height / 2);
        float avgX = (V1.Position.x + V2.Position.x) / 2;
        float avgY = (V1.Position.y + V2.Position.y) / 2;
        float rotation = atan2(V1.Position.y - V2.Position.y, V1.Position.x - V2.Position.x);
        // converts radians to degrees
        rotation = rotation * 180 / 3.14159;
        rect.setPosition(avgX+200, avgY+200);
        rect.setRotation(rotation + 90);
        rect.setFillColor(sf::Color::Blue);
        window.draw(rect);

        // draws a line showing the collision vector
        V1.Position = sf::Vector2f(0, 0);
        V2.Position = physics::CollisionInfo.Normal * physics::CollisionInfo.Depth;
        V2.Position.x *= 5;
        V2.Position.y *= 5;
        rect.setSize(sf::Vector2f(2, distanceBetweenVertices(V1.Position, V2.Position)));
        rect.setOrigin(rect.getGlobalBounds().width / 2, rect.getGlobalBounds().height / 2);
         avgX = (V1.Position.x + V2.Position.x) / 2;
         avgY = (V1.Position.y + V2.Position.y) / 2;
         rotation = atan2(V1.Position.y - V2.Position.y, V1.Position.x - V2.Position.x);
        // converts radians to degrees
        rotation = rotation * 180 / 3.14159;
        rect.setPosition(avgX + 200, avgY + 200);
        rect.setRotation(rotation + 90);
        rect.setFillColor(sf::Color::Red);
        window.draw(rect);
    }

    // complex polygons
    for (unsigned int i = 0; i < ComplexPolygons.size(); i++) {
        ComplexPolygon &polygon = ComplexPolygons[i];
        for (unsigned int t = 0; t < polygon.triangles.size(); t++) {
            vector<Edge> T = polygon.triangles[t];
            for (unsigned int e = 0; e < T.size(); e++) {
                sf::CircleShape circle(4);
                circle.setPosition(T[e].v1->Position);
                circle.setOrigin(4, 4);
                circle.setFillColor(sf::Color::Green);
                window.draw(circle);

                circle.setPosition(T[e].v2->Position);
                window.draw(circle);

                if (T[e].Collision) {
                    Vertex V1 = *T[e].v1;
                    Vertex V2 = *T[e].v2;
                    sf::RectangleShape rect(sf::Vector2f(10 - (t * 2), distanceBetweenVertices(V1.Position, V2.Position)));
                    rect.setOrigin(rect.getGlobalBounds().width / 2, rect.getGlobalBounds().height / 2);
                    float avgX = (V1.Position.x + V2.Position.x) / 2;
                    float avgY = (V1.Position.y + V2.Position.y) / 2;
                    float rotation = atan2(V1.Position.y - V2.Position.y, V1.Position.x - V2.Position.x);
                    // converts radians to degrees
                    rotation = rotation * 180 / 3.14159;
                    rect.setPosition(avgX, avgY);
                    rect.setRotation(rotation + 90);
                    rect.setFillColor(sf::Color(t * 40 + 120, 120, 120));
                    window.draw(rect);
                }
            }
        }
    }

    // circles
    for (unsigned int c = 0; c < physics::circles.size(); c++) {
        physics::Circle& circle = physics::circles[c];
        sf::CircleShape circleShape;
        circleShape.setOrigin(sf::Vector2f(circle.radius, circle.radius));
        circleShape.setPosition(circle.Pos);
        circleShape.setRadius(circle.radius);
        circleShape.setOutlineThickness(2);
        circleShape.setOutlineColor(sf::Color::Green);
        circleShape.setFillColor(sf::Color::Transparent);
        window.draw(circleShape);
    }
}

void RenderNormals(sf::RenderWindow &window) {
    using namespace physics;
    for (unsigned int i = 0; i < SimplePolygons.size(); i++) {
        SimplePolygon& p = SimplePolygons[i];
        for (unsigned int j = 0; j < p.edges.size(); j++) {
            Edge& E = p.edges[j];
            sf::Vector2f temp = E.v1->Position - E.v2->Position;
            sf::Vector2f Normal = sf::Vector2f(temp.y, temp.x);
            Normal = sf::Vector2f(E.v1->Position.y - E.v2->Position.y, E.v2->Position.x - E.v1->Position.x);
            Normal = normalize(Normal);
            Normal.x *= 15;
            Normal.y *= 15;

            // first circle
            sf::CircleShape circle;
            circle.setFillColor(sf::Color::Yellow);
            circle.setRadius(3);
            circle.setOrigin(sf::Vector2f(circle.getRadius(), circle.getRadius()));
            float lineX = (E.v1->Position.x + E.v2->Position.x) / 2;
            float lineY = (E.v1->Position.y + E.v2->Position.y) / 2;

            circle.setPosition(sf::Vector2f(lineX,lineY) + Normal);
            window.draw(circle);

            // second circle
            circle.setPosition(sf::Vector2f(lineX, lineY));
            window.draw(circle);

            sf::Vector2f p1 = sf::Vector2f(lineX, lineY);
            sf::Vector2f p2 = sf::Vector2f(lineX, lineY) + Normal;
            float rotation = atan2(p1.y - p2.y, p1.x - p2.x);
            // converts radians to degrees
            rotation = rotation * 180 / 3.14159;
            sf::RectangleShape rect;
            rect.setPosition(Normal.x/2+lineX, Normal.y/2+lineY);
            rect.setSize(sf::Vector2f(2, distanceBetweenVertices(p1, p2)));
            rect.setOrigin(sf::Vector2f(rect.getGlobalBounds().width / 2, rect.getGlobalBounds().height / 2));
            rect.setRotation(rotation + 90);
            rect.setFillColor(sf::Color::Yellow);
            window.draw(rect);
        }
    }
}

#include <chrono>
#include <thread>


int main() {

    using namespace std;

    if (!font.loadFromFile("LicensePlate.ttf"))
    {
        // error...
    }
   
    //physics::MakeComplexPolygon({sf::Vector2f(70, 0),sf::Vector2f(70,60), sf::Vector2f(120,60), sf::Vector2f(120, 0)});
    //physics::MakeComplexPolygon({ sf::Vector2f(170, 0),sf::Vector2f(170,60), sf::Vector2f(220,60), sf::Vector2f(220, 0) });

    physics::MakePolygon({sf::Vector2f(70,70),sf::Vector2f(70,200), sf::Vector2f(120,200), sf::Vector2f(120, 70)});
    physics::MakePolygon({ sf::Vector2f(170,70),sf::Vector2f(170,200), sf::Vector2f(220,200), sf::Vector2f(220, 70) });

    //physics::MakeComplexPolygon({ sf::Vector2f(70,70), sf::Vector2f(95, 120),sf::Vector2f(95, 175), sf::Vector2f(70,200), sf::Vector2f(120,200), sf::Vector2f(120, 70)});
    //physics::MakeCircle(sf::Vector2f(230, 150), 25);
    //physics::MakeCircle(sf::Vector2f(300, 150), 10);

    /*physics::makePhysicsConstraint(physics::ComplexPolygons[0].vertices[1], physics::ComplexPolygons[1].vertices[0], physics::ComplexPolygons[1].triangles[0][1],
        1, .05, true);*/

    sf::RenderWindow window(sf::VideoMode(800, 800), "Physics Test");

    sf::Clock clock;
    clock.restart();
    while (window.isOpen()) {
        sf::Event event;
        window.clear();
        sf::Vector2i mousePos = sf::Mouse::getPosition(window);
        sf::Vector2f mousePosf = sf::Vector2f(mousePos.x, mousePos.y);
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
            if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
            {
                physics::dragShapes(mousePosf);
            }
            else {
                physics::debuffLeft = false;
            }
        }
        float DELTA_TIME = clock.restart().asSeconds();
        physics::update(DELTA_TIME);
        physics::contrainShapes(physics::SimplePolygons, window);
        RenderNormals(window);
        renderVertices(window);
        window.display();
    }
    return 0;
}
