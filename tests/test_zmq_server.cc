// 
// Hello World server in C++
// Binds REP socket to tcp://*:5555
// Expects "Hello" from client, replies with "World"
//
#include <zmq.hpp>
#include <string>
#include <iostream>
#include <unistd.h>

int main()
{
    using namespace std;

    // Prepare our context and socket
    zmq::context_t context(1);
    zmq::socket_t socket(context, ZMQ_REP);
    socket.bind("tcp://*:5555");
    cout << "Starting Hello World server on port 5555..." << endl;

    while (true)
    {
        zmq::message_t request;

        // Wait for next request from client
        socket.recv(&request);
        string hello = (char *)request.data();
        cout << "Received " << hello << endl;

        // Do some 'work'
        sleep(1);

        // Send reply back to client
        zmq::message_t reply(5);
        memcpy((void *)reply.data(), "World", 5);
        socket.send(reply);
    }

    return 0;
}
