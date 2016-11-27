#ifndef EVENTFILTERMANAGER_H
#define EVENTFILTERMANAGER_H

#include <QWidget>
#include <QMap>
#include <qDebug>

class EventFilterManager : public QObject
{
	Q_OBJECT

public:
	EventFilterManager(QObject* parent = 0);
	~EventFilterManager();

	void addFilterWidget(const QString name, QObject* filter);
	QObject* getFilterWidget(const QString name);
	QObject* removeFilterWidget(const QString name);
	QString const getCurrentFilterName() const;
	void setCurrentFilterName(const QString name);
	QList<QString> getFilterNames();
	void clear();

protected:
	bool eventFilter(QObject* target, QEvent* event);

private:
	QMap<QString, QObject*> filterMap;
	QString currentFilter;
};

#endif // EVENTFILTERMANAGER_H

