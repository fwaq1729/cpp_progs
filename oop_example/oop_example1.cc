#include <string>
#include <iostream>
#include <vector>
using namespace std;

const constexpr static int nvals = 10;
const constexpr static int nvals1 = 20;

class Base_Employee {
  virtual void AskForPromotion() = 0;
};

class Employee : Base_Employee { // Using inheritance, create a derived class
private:
  string Company;
  size_t Age;
  // Keep same order for sizing for Data Integrity
  vector<size_t> arr_things_;
  vector<string> arr_things1_;  

protected:
  string Name; // to give access in derived class

public:
 // Define setters:
  void setName(const string name) {this->Name = name;} // Using this to point to private variables
  void setCompany(const string company) { this->Company = company;}
  void setAge(const size_t age) { this->Age = age; }
  // Define getters:
  string getName () { return Name; }  
  string getCompany() { return Company; }  
  size_t getAge() { return Age; }

  void Introduce_Yourself() {
    cout << "Hi, my name is " << Name << ", I am " << Age << " old. I work in " << Company << endl;
  }

  // tip1: Using const on input variables for safety
  // tip2: Passing values using AAA(BBB) mode
  // tip3: this is a non-default constructor that does nothing (only passes variables)
  Employee(const string name, const string company, const size_t age)
    : Name(name), Company(company), Age(age),
    arr_things_(nvals),  // Keep same order for sizing for Data Integrity
    arr_things1_(nvals)
  {
    cout << "arrays were created, size(arr_things_) = " << arr_things_.size() << endl;
    cout << "arrays were created, size(arr_things1_) = " << arr_things1_.size() << endl;
  }

  void AskForPromotion() override { // Override base method
    if (Age > 35) {
      cout << Name << " got promoted!" << endl;
    } else {
      cout << Name << " does not deseve promotion ..." << endl;
    }
  }

  virtual void Task() {
    cout << Name << " is checking email, task backlog, performing tasks ..." << endl;
  }
};

class Developer : public Employee { // Using inheritance, create a derived class
public:
  string FavProgrammingLanguage;
  Developer(const string name, const string company, const size_t age, const string favProgrammingLanguage)
    : Employee(name, company, age), FavProgrammingLanguage(favProgrammingLanguage)
  {
  }
  void Wrote_Code() {
    cout << Name << " developed code using " << FavProgrammingLanguage << endl;
  }
  void Task() override { // Override base method
    cout << Name << " is writing " << FavProgrammingLanguage << " code " << endl;
  }
};

class Teacher : public Employee { // Using inheritance, create a derived class
public:
  string Subject;
  void PrepareLesson() {
    cout << Name << " is preparing " << Subject << " lesson" << endl;
  }
  Teacher(const string name, const string company, const size_t age, string subject)
    : Employee(name, company, age), Subject(subject)
  {
  }
  void Task() override { // Override base method
    cout << Name << " is teaching " << Subject << " class " << endl;
  }
};

int main()
{
  Developer d = Developer("Roger", "QuantumComputing", 30, "C++");
  cout << "check name1: name is " << d.getName() << endl;
  d.setName("Rose");
  cout << "check name2: name is " << d.getName() << endl;
  Teacher t = Teacher("Luigi", "Montecito High School", 45, "Literature");
  t.PrepareLesson();
  t.AskForPromotion();
  d.Wrote_Code();
  d.AskForPromotion(); // Access method from inherited class Employee
  d.Task();
  t.Task();
  // Using polymorphism, parent class reference is used to refere to a child class object
  Employee *p1 = &d;
  Employee *p2 = &t;
  p1->Task();
  p2->Task();
  return 0;  
}
