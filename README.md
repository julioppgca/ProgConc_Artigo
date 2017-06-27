# Source files of the paper "Communication between Processors in Multiplatform Systems".

This paper presents a practical implementation of the distributed processing concepts applied to the Single Board Computer Udoo Dual Lite. For this purpose, a data acquisition system, based on the Direct Memory Access Controller, is implemented in an ARM Cortex-M3 microcontroller. The collected data is sent continuously, through an UART interface, to a second microcontroller that runs Linux Ubuntu operating system. A software written in C uses the pthread library for synchronization between read and write operations of the received data. A web server provides an HTML page for querying and viewing the collected data signal.

## Additional steps:

Compile main.c with:

  >> gcc main.c -o serial -lpthread -lm
  
Copy the compiled file, "serial", to the "web" folder.

Install "socket.io" and "chart.js"

  >> npm install socket.io --save
  >> npm install chart.js --save
  
Run the server with the command:

  >> node web.js
  
Connecet DAC1 to AN0, check arduino pinout diagram.

Verify the IP of the board:

  >> ifconfig eth0 
  
At the browser insert the IP address with the port number 8081, eg.:

  >> http://192.168.2.7:8081
