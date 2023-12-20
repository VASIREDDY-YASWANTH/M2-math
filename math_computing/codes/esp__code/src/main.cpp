#include <Arduino.h>
#ifdef ESP32
  #include <WiFi.h>
  #include <AsyncTCP.h>
#else
  #include <ESP8266WiFi.h>
  #include <ESPAsyncTCP.h>
#endif
#include <ESPAsyncWebServer.h>
#include"matfun.h"
#include"libr.h"

AsyncWebServer server(80);

const char* ssid = "userv";
const char* password = "wordabcd";

const char* input_parameter00 = "input00";
const char* input_parameter01 = "input01";
const char* input_parameter10 = "input10";
const char* input_parameter11 = "input11";
const char* input_parameter20 = "input20";
const char* matrix1[2]={input_parameter00,input_parameter01};     // matrix for vertex A
const char* matrix2[2]={input_parameter10,input_parameter11};     // matrix for incentre B
const char* matrix3[1]={input_parameter20};

const char index_html[] PROGMEM = R"rawliteral(
<!DOCTYPE HTML><html><head>
  <title>CONTACT POINTS</title>
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <style>
    html {font-family: Times New Roman; display: inline-block; text-align: center;}
    h2 {font-size: 2.0rem; color: blue;}
  </style> 
  </head><body>
  <h2>TO FIND CONTACT POINTS USING EIGEN APPROACH</h2> 
  <p>Enter the values of points A,B,C
  <form action="/get">
    Enter the values of Point A: <input type="number" name="input00"><br><input type="number" name="input01"><br>
    Enter the values of Point B: <input type="number" name="input10"><br><input type="number" name="input11"><br>
    Enter the value of radius r: <input type="number" name="input20"><br>
    <input type="submit" value="Submit">
    

  </form><br>
</body></html>)rawliteral";

void notFound(AsyncWebServerRequest *request) {
  request->send(404, "text/plain", "Not found");
}

void setup() {
  Serial.begin(115200);
  WiFi.mode(WIFI_STA);
  WiFi.begin(ssid, password);
  if (WiFi.waitForConnectResult() != WL_CONNECTED) {
    Serial.println("Connecting...");
    return;
  }
  Serial.println();
  Serial.print("IP Address: ");
  Serial.println(WiFi.localIP());

  server.on("/", HTTP_GET, [](AsyncWebServerRequest *request){
    request->send_P(200, "text/html", index_html);
  });

  server.on("/get", HTTP_GET, [] (AsyncWebServerRequest *request) {
   double **A,**B,**C;
    A=load_ser(request,matrix1,2);
    B=load_ser(request,matrix2,2);
    C=load_ser(request,matrix3,1);
    
    double **h=createMat(2,1);     //vertex    
    double **u=createMat(2,1);	   //incenter
    double **V=Mateye(2);	   //2x2 identity matrix
    h[0][0]=A[0][0];    h[1][0]=A[1][0];
    u[0][0]=B[0][0];    u[1][0]=B[1][0];
    double r=C[0][0];  		   //inradius
    double **f=createMat(1,1); f[0][0]=pow(Matnorm(u,2),2)-r*r;

	double **cp=contactPoints(V,h,u,f); 	//function for finding contact points 

  request->send(200, "text/html", "Contact point Matrix: <br>"+ String(cp[0][0],6) + ","+ String(cp[1][0],6)+"<br>"+ String(cp[0][1],6)+","+ String(cp[1][1],6)+  ")<br><a href=\"/\">Return to Home Page</a>");


  });
  server.onNotFound(notFound);
  server.begin();
}

void loop() {
  
}
