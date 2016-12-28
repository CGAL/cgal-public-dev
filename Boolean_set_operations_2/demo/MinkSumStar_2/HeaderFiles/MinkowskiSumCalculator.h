#ifndef CGAL_MINKOWSKISUMCALCULATOR_H
#define CGAL_MINKOWSKISUMCALCULATOR_H

#include <list>
#include <map>
//#include <memory>

#include <QThread>

#include "Typedefs.h"

class MinkowskiSumResult
{
public:
  MinkowskiSumResult();
  MinkowskiSumResult(const StdPolygonList& sum);

  const StdPolygonList& getSum() const { return sum_; }
  const std::list<Polygon>& getDecomposedSum() const { return decomposedSum_; }
  double getArea() const;
  int getSize() const;
  void print() const;

private:
  StdPolygonList sum_;
  std::list<Polygon> decomposedSum_;
};

class MinkowskiSumCalculator : public QObject
{
  Q_OBJECT

public:
  MinkowskiSumCalculator();
  ~MinkowskiSumCalculator();

  int getKValue();

  void setInput(const StdPolygonList& input);
  std::shared_ptr<const MinkowskiSumResult> getSum(int k);

  static void decomposePolygon(const PolygonWithHoles& polygon, std::list<Polygon>& result);

private:
  PolygonWithHoles scalePolygon(const PolygonWithHoles& polygon, Kernel::FT scale);
  PolygonWithHoles removeCollinearPoints(const PolygonWithHoles& polygon);
  Polygon removeCollinearPoints(const Polygon& polygon);

  std::map<int, std::shared_ptr<MinkowskiSumResult>> cache_;
  StdPolygonList input_;
  int kValue_;

public slots:
  void nextSum();
  void toSum(int k);
  void slotSetInput(const StdPolygonList& input);

signals:
  void polygonReady(StdPolygonList* polygon, int k);
};

class MinkowSkiSumStarController : public QObject
{
  Q_OBJECT
    QThread m_workerThread;
  MinkowskiSumCalculator* m_minkSum;
  bool m_running;
  bool m_inputIsSet;
  int m_maxK;

public:
  MinkowSkiSumStarController();
  ~MinkowSkiSumStarController();

  void setInput(const StdPolygonList& input);
  void runCalculation(int maxK = -1);

public slots:
  void handleResults(StdPolygonList* polygon, int k);

signals:
  void operate();
  void signalSetInput(const StdPolygonList& input);
  void requestNewPolygonCreation(PolygonWithHoles* polygon, QString name, int layer);
};

#endif // CGAL_MINKOWSKISUMCALCULATOR_H
