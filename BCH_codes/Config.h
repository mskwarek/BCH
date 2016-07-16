class Config{
 private:
  std::map<std::string, std::string> databaseInfo;
 public:
  void load(const std::string&);
  std::map<std::string, std::string> getDatabaseInfo();
};
