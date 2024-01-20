const int speakerPins[] = {10, 11, 8, 9, 6, 7, 4, 5, 13, 12, A4, A5, A1, A0, A2, A3};
const int numSpeakers = 16;

void setup() {
  Serial.begin(115200); // Start serial communication
  Serial.setTimeout(100);
  for (int i = 0; i < numSpeakers; i++) {
    pinMode(speakerPins[i], OUTPUT); // Set each speaker pin as output
    digitalWrite(speakerPins[i], LOW); // Turn all speakers off initially
  }
}

void loop() {

  if (Serial.available()) {
    String input = Serial.readStringUntil('>');
    if (input.startsWith("<")) {
      input = input.substring(1, input.length()); // Remove start and end markers
        int colonIndex = input.indexOf(":");
        if (colonIndex > 0) {
          int state = input.substring(0, colonIndex).toInt(); //turn on or off
          input = input.substring(colonIndex+1); //first number after the colon
          
          while (input.length() > 0) {
            int commaIndex = input.indexOf(",");
            if (commaIndex < 0) {
              commaIndex = input.length();
            }
            int speakerIndex = input.substring(0, commaIndex).toInt() - 1; // Speaker indexes start from 1, so subtract 1
            if (speakerIndex >= 0 && speakerIndex <= numSpeakers) { 
              digitalWrite(speakerPins[speakerIndex], state); 
            }
            input = input.substring(commaIndex+1);
        }
      }
    }
  }
}
