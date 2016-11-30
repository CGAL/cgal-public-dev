#include "EventFilterManager.h"

EventFilterManager::EventFilterManager(QObject* parent) :
	QObject(parent),
	filterMap()
{
	QString currentFilter("");
}

EventFilterManager::~EventFilterManager()
{
}

void EventFilterManager::addFilterWidget(QString name, QObject* filter)
{
	filterMap[name] = filter;
}

QObject* EventFilterManager::getFilterWidget(QString name)
{
	return filterMap[name];
}

QObject* EventFilterManager::removeFilterWidget(QString name)
{
	QObject* filter = filterMap[name];
	filterMap.remove(name);
	if (currentFilter == name) {
		currentFilter = "";
	}
	return filter;
}

QString const EventFilterManager::getCurrentFilterName() const
{
	return currentFilter;
}

void EventFilterManager::setCurrentFilterName(QString name)
{
	if (filterMap.contains(name) || name == "") {
		currentFilter = name;
	}
	else {
		qDebug() << Q_FUNC_INFO << "The name \"" + name + "\" is not a name of any filter and is not an empty name";
	}
}

QList<QString> EventFilterManager::getFilterNames()
{
	return filterMap.uniqueKeys();
}

void EventFilterManager::clear()
{

	currentFilter = "";
	filterMap.clear();
}

bool EventFilterManager::eventFilter(QObject* target, QEvent* event)
{
	if (filterMap.contains(currentFilter)) {
		return filterMap[currentFilter]->eventFilter(target, event);
	}
	return false;
}
