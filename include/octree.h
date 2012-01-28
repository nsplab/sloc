#ifndef SLOC_OCTREE_H
#define SLOC_OCTREE_H

#include <vector>

namespace sloc {
namespace octree {
// ----------------------------------------------------------------------------

class Point
{
public:
    long id;
    double x, y, z;
    unsigned int code;

public:
    Point(long id, double x, double y, double z)
        : id(id), x(x), y(y), z(z)
    {
    }

    Point(const Point& pt)
    {
        id = pt.id; 
        x = pt.x; y = pt.y; z = pt.z;
        code = pt.code;
    }
};

struct Bounds
{
    Point center;
    double size;
    //double maxx,maxy,maxz;
    //double minx,miny,minz;

    Bounds(Point center, double size)
        : center(center), size(size)
    {
        //calc_min_max();
    }

    Bounds(const Bounds& b)
    {
        center = b.center;
        size = b.size;
        //maxx = b.maxx; maxy = b.maxy; maxz = b.maxz;
        //minx = b.minx; miny = b.miny; minz = b.minz;
    }

    bool contains(const Point& p)
    {
        double maxx,maxy,maxz;
        double minx,miny,minz;
        maxx = center.x + size/2;
        maxy = center.y + size/2;
        maxz = center.z + size/2;
        minx = center.x - size/2;
        miny = center.y - size/2;
        minz = center.z - size/2;
        return ((minx <= p.x) && (p.x <= maxx)) &&
               ((miny <= p.y) && (p.y <= maxy)) &&
               ((minz <= p.z) && (p.z <= maxz));
    }

    void set_child_bounds(int n, Bounds& cb)
    {
        //
        // Calculate bits of n
        //
        //   i = n % 2;
        //   j = (n/2) % 2;
        //   k = (n/4) % 2;
        //

        int i = n & 1;
        int j = (n >> 1) & 1;
        int k = (n >> 2) & 1;

        //
        // Select an octant based on the bit pattern of n
        //
        cb.center.x = center.x + (2*i-1) * size/4;
        cb.center.y = center.y + (2*j-1) * size/4;
        cb.center.z = center.z + (2*k-1) * size/4;

        //cb.calc_min_max();
    }
};

struct Node
{
    Bounds bounds;
    Node* children[8];
    std::vector<Point> points;

    Node()
    {
        for (int i = 0; i < 8; i++)
            children[i] = 0;
    }

    ~Node()
    {
        clear();
    }

    void clear()
    {
        for (int i = 0; i < 8; i++)
        {
            if (children[i] != 0)
            {
                children[i]->clear();
                delete children[i];
                children[i] = 0;
            }
        }
    }

    void add(long id, double x, double y, double z, int level)
    {
        // determine to which of the children
        // our point falls into

        Node* child = 0;


        if ((level >= maxlevel) || (child == 0))
        {
            // either we've reached the max recursion level,
            // or the point belongs to more than one child
            // (straddling on boundary).

            Point p(id,x,y,z);

            // XXX: check that point is not already on our list

            points.push_back(p);
        }
        else (child != 0)
        {
            // our point belongs to only one child.
            // recursively add it to that node.
            child->add(level+1, x, y, z, id);
        }
    }

    bool find(long& id, double x, double y, double z)
    {
        int i;

        // search points
        for (i = 0; i < points.size(); i++)
        {
            // search points
        }

        // search children

        return false;
    }

    int code(double x, double y, double z)
    {
    }

    size_t bytes()
    {
        size_t b = sizeof(Node);
        b += points.size() * sizeof(Point);
        for (int i = 0; i < 8; i++)
        {
            if (children[i] != 0)
                b += children[i]->bytes();
        }
        return b;
    }
};


class Octree
{
public:
    Octree(int levels)
    {
        maxlevel = levels;
    }

    ~Octree()
    {
    }

    void add(long id, double x, double y, double z)
    {
        root.add(id,x,y,z);
    }

    bool find(long& id, double x, double y, double z)
    {
        return root.find(id,x,y,z);
    }

    size_t bytes() const
    {
        return root.bytes();
    }

public:
    int maxlevel;
    Node root;
};

// ----------------------------------------------------------------------------
}} // namespace sloc::octree

#endif
