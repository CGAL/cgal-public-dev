#ifndef EVENTFILTERMANAGERGROUP_H
#define EVENTFILTERMANAGERGROUP_H

#include <QWidget>
#include <QMap>
#include <qDebug>

#include "EventFilterManager.h"

class EventFilterManagerGroup
{

public:
	EventFilterManagerGroup();
	~EventFilterManagerGroup();

	void addObjectToWatch(const QString name, QObject* object);
	void removeObjectFromWatch(const QString name);
	void addFilterWidget(const QString targetName, const QString name, QObject* filter);
	QObject* removeFilterWidget(const QString name);
	QObject* getFilterWidget(const QString name);
	QString const getCurrentFilterName() const;
	void setCurrentFilterName(const QString name);
	QList<QString>* getFilterNames();
	void clear();

private:
	QMap<QString, EventFilterManager*> objectsMap;
	QMap<QString, QString> filtersMap;
	QString currentFilter;
	QString currentWatchedObject;

};

#endif // EVENTFILTERMANAGERGROUP_H
