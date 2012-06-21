//
// Hello World client in C++
// Connects REQ socket to tcp://localhost:5555
// Sends "Hello" to server, expects "World" back
//
#include <zmq.hpp>
#include <string>
#include <iostream>

int main()
{
    using namespace std;

    // Prepare our context and socket
    zmq::context_t context(1);
    zmq::socket_t socket(context, ZMQ_REQ);

    cout << "Connecting to hello world server..." << endl;
    socket.connect("tcp://localhost:5555");

    // Do 10 requests, waiting each time for a response
    for (int n = 0; n < 10; n++)
    {
        zmq::message_t request(6);
        memcpy((void *)request.data(), "Hello", 5);
        cout << "Sending Hello " << n << "... " << std::flush;
        socket.send(request);

        // Get the reply
        zmq::message_t reply;
        socket.recv(&reply);
        string world = (char *)reply.data();
        cout << "-> Received " << world << " " << n << endl;
    }

    return 0;
}
