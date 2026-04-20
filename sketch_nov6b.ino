#include <Servo.h>

Servo esc_1;

// ====== USER SETTINGS ======
int powerLevel = 100;             // -100..100
unsigned long runTime  = 2800;    // ms for each direction
unsigned long idleTime = 500;     // ms stop time between directions

// ====== Encoder pins ======
const int encoderA = 2;  // D2 (interrupt)

// ====== Encoder counters / timing ======
volatile long encoderCount = 0;     // counts from Channel A only
unsigned long lastReportMs = 0;
const unsigned long REPORT_PERIOD_MS = 100;  // ms
const double COUNTS_PER_REV = 103.8;         // single-channel pulses per rev

// ====== Motor timing (state machine) ======
enum MotorState { FORWARD, STOP1, BACKWARD, STOP2 };
MotorState motorState = FORWARD;
unsigned long stateStartMs = 0;

void setup() {
  // ESC
  esc_1.attach(9);
  delay(2000);  // arm ESC

  // Encoder setup
  pinMode(encoderA, INPUT_PULLUP);
  attachInterrupt(digitalPinToInterrupt(encoderA), isrA, RISING);  // count rising edges

  // Serial (CSV only)
  Serial.begin(9600);
  while (!Serial) { ; }

  // Print CSV header once (no other text)
  Serial.println("RPM,TotalCount");

  stateStartMs = millis();  // start timing
}

void loop() {
  unsigned long now = millis();

  // ----------------- Compute RPM & print CSV -----------------
  if (now - lastReportMs >= REPORT_PERIOD_MS) {
    lastReportMs = now;

    // snapshot encoder counts safely
    noInterrupts();
    long countNow = encoderCount;
    interrupts();

    static long lastCount = 0;
    long delta = countNow - lastCount;
    lastCount = countNow;

    // convert to RPM
    double revs = delta / COUNTS_PER_REV;
    double rpm = revs * (60000.0 / REPORT_PERIOD_MS);  // 60s/min × 1000ms/s

    // CSV output: RPM,TotalEncoderCount
    Serial.print(rpm, 2);
    Serial.print(",");
    Serial.println(countNow);
  }

  // ----------------- Motor non-blocking state machine -----------------
  switch (motorState) {
    case FORWARD:
      set_esc_power(powerLevel);
      if (now - stateStartMs >= runTime) {
        motorState = STOP1;
        stateStartMs = now;
      }
      break;

    case STOP1:
      set_esc_power(0);
      if (now - stateStartMs >= idleTime) {
        motorState = BACKWARD;
        stateStartMs = now;
      }
      break;

    case BACKWARD:
      set_esc_power(-powerLevel);
      if (now - stateStartMs >= runTime) {
        motorState = STOP2;
        stateStartMs = now;
      }
      break;

    case STOP2:
      set_esc_power(0);
      if (now - stateStartMs >= idleTime) {
        motorState = FORWARD;
        stateStartMs = now;
      }
      break;
  }
}

// ----------------- Helpers -----------------
void set_esc_power(int power) {
  power = constrain(power, -100, 100);
  const int signal_min = 1050;
  const int signal_max = 1950;
  int signal_output = map(power, -100, 100, signal_min, signal_max);
  esc_1.writeMicroseconds(signal_output);
}

// ----------------- Encoder ISR -----------------
void isrA() {
  encoderCount++;  // increment on every rising edge of channel A
}